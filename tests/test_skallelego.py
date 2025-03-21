import allel
import skallelego
import os
import pytest

# Read input VCF with scikit-allel

test_dir = os.path.join(os.path.dirname(__file__))
test_vcf = os.path.join(test_dir, "data", "Pundamilia.RAD.vcf.gz")
test_output = os.path.join(test_dir, "data", "output.vcf.gz")


fields = [
    'samples', 
    'calldata/GT', 
    'variants/ALT', 
    'variants/CHROM', 
    'variants/ID', 
    'variants/POS', 
    'variants/QUAL', 
    'variants/REF',
    'variants/is_snp'
]

@pytest.fixture()
def test_variant_filter():
    callset = allel.read_vcf(test_vcf, fields=fields)
    is_snp = callset['variants/is_snp'][:] #np.ndarray
    ref = callset['variants/REF'][:].astype('U13') #np.ndarray
    alt = callset['variants/ALT'][:].astype('U13') #list
    qual = callset['variants/QUAL'][:]
    gt = allel.GenotypeArray(callset['calldata/GT'][:])
    samples = callset['samples'][:].astype('U13')
    pos = callset['variants/POS'][:]
    chrom = callset['variants/CHROM'][:].astype('U13')

    chrom_flt, pos_flt, ref_flt, alt_flt, qual_flt, gt_flt = skallelego.variant_filter(
        gt,
        chrom,
        pos,
        ref,
        alt,
        qual,
        is_snp,
        core_variants=True,
        ld_threshold=0.2
    )

    assert gt_flt.shape[0] == 41

    return (chrom_flt, pos_flt, ref_flt, alt_flt, qual_flt, gt_flt, samples)

    

def test_write_vcf(test_variant_filter):
    chrom_flt, pos_flt, ref_flt, alt_flt, qual_flt, gt_flt, samples = test_variant_filter
    success = skallelego.write_vcf(
        samples,
        chrom_flt,
        pos_flt,
        ref_flt,
        alt_flt,
        qual_flt,
        gt_flt,
        os.path.join(test_dir, "output.vcf.gz")
    )
    assert success == True
    assert os.path.exists(os.path.join(test_dir, "output.vcf.gz")) == True


def test_clean_test():
    test_output = os.path.join(test_dir, "output.vcf.gz")
    os.remove(test_output)
    assert os.path.exists(test_output) == False
