# skallelego

**Scikit-allel lego**: a few utility functions to complement [scikit-allel](https://github.com/cggh/scikit-allel)

Write scikit-allel objects to disk as VCF(.gz)

## Installation

```bash
pip install git+https://github.com/maxibor/skallelego.git
```

## Quick start

```python
import scikit-allel
import skallelego

# Read input VCF with scikit-allel
callset = allel.read_vcf('example.vcf')
is_snp = callset['variants/is_snp'][:] #np.ndarray
ref = callset['variants/REF'][:].astype('U13') #np.ndarray
alt = callset['variants/ALT'][:].astype('U13') #list
qual = callset['variants/QUAL'][:] if "QUAL" in callset['variants'] else None  #np.ndarray
gt = allel.GenotypeArray(callset['calldata/GT'][:])
samples = callset['samples'][:].astype('U13')
pos = callset['variants/POS'][:]
chrom = callset['variants/CHROM'][:].astype('U13')

# Do some filtering/variants operations
# ...

# Write objects to disk
skallelego.vcf2disk(
    samples,
    chrom,
    pos,
    ref,
    alt,
    qual,
    gt,
    "output.vcf.gz"
)
```

## Documentation

```python
>>>help(skallelego.vcf2disk)
Help on function vcf2disk in module skallelego:

vcf2disk(samples, chrom, pos, ref, alt, qual, gt, output)
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

variant_filter(gt_in, chrom, pos, ref, alt, qual)
    Filter variants to only keep non singleton bi-allelic variants, free from LD
    Args:
        gt_in(allel.model.ndarray.GenotypeArray): scikit-allel genotype array of shape (V,S,P), P being the ploidy
        chrom(np.array): array of chromosome names of shape (V,)
        pos(np.array): array of positions of shape (V,)
        ref(np.array): array of reference alleles of shape (V,)
        alt(np.array): array of alternate alleles of shape (V, N) (N being the max number of alternate allels)
        qual(np.array): array of variant qualities of shape (V,)
    
    Returns:
        np.array: filtered chromosome array, of shape (Vf, )
        np.array: filtered position array, of shape (Vf, )
        np.array: filtered reference allele array, of shape (Vf, )
        np.array: filtered alternate allele array, of shape (Vf, 2)
        np.array: filtered variant quality array of shape (Vf, )
        allel.model.ndarray.GenotypeArray: filtered genotype array of shape (Vf, S, 2)
```
