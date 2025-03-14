import vcfpy
import allel
import numpy as np
from tqdm import tqdm

def ld_prune(gn, size=500, step=200, threshold=.1, n_iter=2):
    """
    Prune variant in Linkage Desequilibrium
    Snippet adapted from http://alimanfoo.github.io/2015/09/28/fast-pca.html

    Args:
        gn(allel.model.ndarray.GenotypeArray): scikit-allel genotype array of shape (V,S,P), P being the ploidy
        size(int): Window size (number of variants).
        step(int): Number of variants to advance to the next window.
        threshold(float): Maximum value of r**2 to include variants.
        n_inter(int): number of iterations

    Returns:
        allel.model.ndarray.GenotypeArray: LD filtered genotype array of shape (Vf, S, P)
        list[np.array]: List (of lenght n_iter) of np.array containing indices of positions kept
    """

    kept_pos = []
    for i in range(n_iter):
        _ = []
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        _.append(list(loc_unlinked))
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
        kept_pos.append(_[0])
    return gn, kept_pos

def variant_filter(gt_in, chrom, pos, ref, alt, qual, is_snp, ld_size=500, ld_step=200, ld_threshold=0.1, ld_iter=2):
    """
    Filter variants to only keep non singleton bi-allelic variants, free from LD
    Args:
        gt_in(allel.model.ndarray.GenotypeArray): scikit-allel genotype array of shape (V,S,P), P being the ploidy
        chrom(np.array): array of chromosome names of shape (V,)
        pos(np.array): array of positions of shape (V,)
        ref(np.array): array of reference alleles of shape (V,)
        alt(np.array): array of alternate alleles of shape (V, N) (N being the max number of alternate allels)
        qual(np.array): array of variant qualities of shape (V,)
        is_snp(np.array): array of binary values for each variant, of shape (V,)
        ld_size(int): Window size (number of variants).
        ld_step(int): Number of variants to advance to the next window.
        ld_threshold(float): Maximum value of r**2 to include variants.
        ld_iter(int): number of iterations
    
    Returns:
        np.array: filtered chromosome array, of shape (Vf, )
        np.array: filtered position array, of shape (Vf, )
        np.array: filtered reference allele array, of shape (Vf, )
        np.array: filtered alternate allele array, of shape (Vf, 2)
        np.array: filtered variant quality array of shape (Vf, )
        allel.model.ndarray.GenotypeArray: filtered genotype array of shape (Vf, S, 2)
    """
    
    gt = gt_in.copy()
    ac = gt.count_alleles()

    # Apply filters
    print(f"{gt.shape[0]} variants - {gt.shape[1]} samples")
    biallelic_mask = ac.is_biallelic()  # Only bi-allelic variants
    non_singleton_mask = ~ac.is_singleton(0) & ~ac.is_singleton(1)  # Exclude singletons
    snp_mask = is_snp  # Only SNPs
    final_mask = snp_mask & biallelic_mask & non_singleton_mask
    print(f"{snp_mask.sum()} SNP are variants")
    print(f"{biallelic_mask.sum()} variants are bi-allelic")
    print(f"{non_singleton_mask.sum()} variants are not singleton")
    nb_variants_remain = final_mask.sum()
    perc_variant_remain = (nb_variants_remain / gt.shape[0])*100
    print(f"Before LD filtering, {round(perc_variant_remain, 2)}% ({nb_variants_remain}) variants remain")



    # Combined mask: SNPs, Bi-allelic, and Non-Singleton
    

    # Filter variants based on the mask
    filtered_chrom = chrom[final_mask]
    filtered_pos = pos[final_mask]
    filtered_ref = ref[final_mask]
    filtered_alt = alt[final_mask]
    filtered_qual = qual[final_mask]
    filtered_gt = gt[final_mask]

    # Perform LD pruning
    pruned_gt, non_ld_mask = ld_prune(filtered_gt.to_n_alt(), size=ld_size, step=ld_step, threshold=ld_threshold, n_iter=ld_iter)

    # Apply LD mask to variants
    for m in non_ld_mask:
        filtered_chrom = filtered_chrom[m]
        filtered_pos = filtered_pos[m]
        filtered_ref = filtered_ref[m]
        filtered_alt = filtered_alt[m][:,:2]
        filtered_qual = filtered_qual[m]
        filtered_gt = filtered_gt[m]

    print(f"After LD filtering {filtered_gt.shape[0]} variants remain")

    return filtered_chrom, filtered_pos, filtered_ref, filtered_alt, filtered_qual, filtered_gt



def write_vcf(samples, chrom, pos, ref, alt, qual, gt, output):
    """
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
    """

    print(f"Writing VCF to disk as {output}")
    
    # Define VCF Header
    s = vcfpy.header.SamplesInfos(samples)
    header = vcfpy.Header()
    header.samples = s

    # Add metadata
    header.add_line(vcfpy.HeaderLine(key="fileformat", value="VCFv4.2"))
    header.add_filter_line({'ID': 'PASS', 'Description': 'All filters passed'})
    header.add_format_line({'ID': 'GT', 'Number': '1', 'Type': 'String', 'Description': 'Genotype'})
    with vcfpy.Writer.from_path(output, header) as writer:
        for i, c in tqdm(enumerate(chrom), total=len(chrom)):
            alt_alleles = [vcfpy.Substitution("SNV", a) for a in alt[i] if a != ""]

            calls = []
            for j, s in enumerate(samples):
                gt_string = "/".join(list(gt[i,j].astype(str)))  
                calls.append(vcfpy.Call(s, {"GT": gt_string}))

            record = vcfpy.Record(
                CHROM=str(c), 
                POS=int(pos[i]),
                ID=".", 
                REF=ref[i],
                ALT=alt_alleles,
                QUAL=float(qual[i]) if qual is not None else None,
                FILTER=["PASS"],  
                INFO={},  
                FORMAT=["GT"],
                calls=calls,
            )

            # Write the record
            writer.write_record(record)