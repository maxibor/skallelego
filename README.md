# skallelego

**Scikit-allel lego**: a few utility functions to complement [scikit-allel](https://github.com/cggh/scikit-allel)

Write scikit-allel objects to disk as VCF(.gz)

## Installation

```bash
pip install git+https://github.com/maxibor/skallelego.git
```

## Quick start

```python
import allel
import skallelego

# Read input VCF with scikit-allel

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

callset = allel.read_vcf('tests/data/Pundamilia.RAD.vcf.gz', fields=fields)
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
    is_snp
)

# Write objects to disk
skallelego.write_vcf(
    samples,
    chrom_flt,
    pos_flt,
    ref_flt,
    alt_flt,
    qual_flt,
    gt_flt,
    "output.vcf"
)
```

## Documentation

```python
>>>help(skallelego.write_vcf)
Help on function write_vcf in module skallelego:

write_vcf(samples, chrom, pos, ref, alt, qual, gt, output)
    Write VCF data to disk. 
    This function currently only supports writing bi-allelic SNVs to VCF files.
    
    Args:
        samples (np.array): List of sample names of shape (S,)
        chrom (np.array): List of chromosome identifiers of shape (V,)
        pos (np.array): List of positions of shape (V,)
        ref (np.array): List of reference alleles of shape(V,)
        alt (np.array): List of alternate alleles of shape (V,P). P is currently set to 2
        qual (np.array): List of variant quality scores of shape(V,)
        gt (scikit-allel genotype array): Genotype matrix of shape(V,S,P)
        output(str): path to output file. Use .gz extension to write bgzf compressed VCF files
    
    Returns:
        None

>>> help(skallelego.variant_filter)
Help on function variant_filter in module skallelego:

variant_filter(gt_in, chrom, pos, ref, alt, qual, is_snp)
    Filter variants to only keep non singleton bi-allelic variants, free from LD
    Args:
        gt_in(allel.model.ndarray.GenotypeArray): scikit-allel genotype array of shape (V,S,P), P being the ploidy
        chrom(np.array): array of chromosome names of shape (V,)
        pos(np.array): array of positions of shape (V,)
        ref(np.array): array of reference alleles of shape (V,)
        alt(np.array): array of alternate alleles of shape (V, N) (N being the max number of alternate allels)
        qual(np.array): array of variant qualities of shape (V,)
        is_snp(np.array): array of binary values for each variant, of shape (V,)

    Returns:
        np.array: filtered chromosome array, of shape (Vf, )
        np.array: filtered position array, of shape (Vf, )
        np.array: filtered reference allele array, of shape (Vf, )
        np.array: filtered alternate allele array, of shape (Vf, 2)
        np.array: filtered variant quality array of shape (Vf, )
        allel.model.ndarray.GenotypeArray: filtered genotype array of shape (Vf, S, 2)
```
