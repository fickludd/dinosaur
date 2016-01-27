import sys
import os
import csv
import numpy
import pymzml
import xml
from pyteomics import mzid


MS1_Precision = 1e-5

def load_feature_table(fn):
    table = []
    with open(fn, 'r') as fh:
        rd = csv.reader(fh, delimiter=',')
        for row in rd:
            if row[0] == 'FEATURE':
                _, rt, mz, _, chg, _, _, _, _, rtl, rtr = row
                table.append([float(mz), int(chg), float(rtl), float(rtr), float(rt)])
    table.sort(key=lambda x: x[3])
    return table


def load_dino_table(fn, minCharge=2):
    table = []
    with open(fn, 'r') as fh:
        rd = csv.reader(fh, delimiter='\t')
        for row in rd:
            if row[0] != 'mz':
                mz, _, chg, rtl, rt, rtr, _, _, _, _, _, _, _, _ = row
                z = int(chg)
                if z >= minCharge:
                    table.append([float(mz), z, float(rtl)*60, float(rtr)*60, float(rt)*60])
    table.sort(key=lambda x: x[3])
    return table


def load_mzid(fn, qval=0.001):
    from pprint import pprint
    psms = []
    specids = [0]
    psmReader = mzid.read(fn)
    for psm in psmReader:
        if psm.has_key('SpectrumIdentificationItem'):
            try:
                specids.append( int(psm['scan number(s)']))
            except KeyError:
                specids.append( int( psm['spectrumID'].split('=')[-1] ))
            else:
                pass

            for match in psm['SpectrumIdentificationItem']:
                if match['MS-GF:QValue'] < qval and match['rank'] == 1 and match['IsotopeError'] == 0 and 2 <= match['chargeState'] <= 4:
                    dm = match['experimentalMassToCharge'] - match['calculatedMassToCharge']
                    dm = dm * 1e6 / match['calculatedMassToCharge']
                    psms.append(dm)
    return numpy.array(psms), max(specids)


def spectra_clone(feature_fn, mzml_fn, dm_offset, max_scan=0, full_iso_width=4.0):
    if feature_fn[-13:] == ".features.csv":
        features = load_dino_table(feature_fn)
    else:
        features = load_feature_table(feature_fn)
    iso_width = full_iso_width / 2.0
    sys.stderr.write("Auto correct precursor m/z offset: %.2f ppm \n" % dm_offset)
    
    if mzml_fn.endswith('.gz'):
        fh = gzip.open(mzml_fn)
    else:
        fh = open(mzml_fn)
    
    outpath = "%s.demix.mgf" % mzml_fn
    sys.stdout = open(outpath, 'wb')

    speciter = pymzml.run.Reader(mzml_fn)
    timescale = 0
    features.sort(key=lambda x: x[0])
    fmz_all = numpy.array([f[0] for f in features])
    try:
        for spec in speciter:
            element = spec.xmlTree.next()
            title = element.get('id')
            idx = int(title.split('scan=')[-1])
            if idx % 1000 == 0 and max_scan > 0:
                sys.stderr.write("DeMix %d MS/MS (~%.1f%%)\n" % (idx, idx * 100.0 / max_scan))

            if not timescale:
                xmltext = xml.etree.ElementTree.tostring(element)
                if xmltext.count(r'unitName="second"'):
                    timescale = 1
                else:
                    timescale = 60

            if spec['ms level'] == 2.0:
                try:
                    rt = float(spec['scan time']) * timescale
                except:
                    continue

                for p in spec['precursors']:
                    pmz = float(p['mz'])
                    try:
                        pz = int(p['charge'])
                    except:
                        pz = 0

                    featured = False
                    peaks = sorted(filter(lambda x: x[1], spec.centroidedPeaks), key=lambda i: i[0])

                    l_idx = fmz_all.searchsorted(pmz - iso_width, side='left')
                    r_idx = fmz_all.searchsorted(pmz + iso_width, side='right')
                    for f in features[l_idx:r_idx]:
                        fmz, fz, frt_left, frt_right, frt = f
                        if frt_left < rt < frt_right:
                            if abs(pmz - fmz) / pmz <= MS1_Precision: 
                                featured = True
                            print 'BEGIN IONS'
                            print 'TITLE=%d[%d:%f:%f]' % (idx, features.index(f), fmz, frt)
                            print 'RTINSECONDS=%f' % rt
                            print 'PEPMASS=%f' % (fmz - fmz * dm_offset * 1e-6)
                            print 'CHARGE=%d+' % fz
                            print 'RAWFILE=%s [%f:%d] diff:%f' % (title, pmz, pz, (fmz - pmz))
                            for a, b in peaks:
                                print a, b
                            print 'END IONS\n'

                    if featured == False and pz > 1:
                        print 'BEGIN IONS'
                        print 'TITLE=%d[-:%f:%f]' % (idx, pmz, rt)
                        print 'RTINSECONDS=%f' % rt
                        print 'PEPMASS=%f' % (pmz - pmz * dm_offset * 1e-6)
                        print 'CHARGE=%d+' % pz
                        print 'RAWFILE=%s' % (title)
                        for a, b in peaks:
                            print a, b
                        print 'END IONS\n'
    except (KeyError, ValueError):
        pass    
    return outpath


if __name__ == '__main__':

    feature_fn = sys.argv[1]    # feature csv table exported from FeatureXML by TOPP. 
    mzml_fn = sys.argv[2]       # centroided MS/MS spectra in mzML, the same file which has been used in the first-pass database search.
    rawpsm_fn = sys.argv[3]     # first-pass database search result: Morpheus .PSMs.tsv file. 
    full_iso_width = float(sys.argv[4]) # the total width of precursor isolation window.
    macc, max_scan = load_mzid(rawpsm_fn)
    # sys.stderr.write("Mean Mass Error (ppm): %.3f SD: %.3f\n" % (macc.mean(), macc.std()))
    # spectra_clone(feature_fn, mzml_fn, macc.mean(), max_scan, full_iso_width)


