/**
 *
 */
package org.theseed.genome.contigs;

import java.util.HashSet;

import org.apache.commons.lang3.StringUtils;

/**
 * This class looks at a position in a contig and accepts or rejects it depending on whether the
 * codon at the current position matches one of several predefined patterns.  The patterns are
 * specified in the constructor.
 *
 * @author Bruce Parrello
 *
 */
public class CodonFilter {

    /** set of codons for which to filter */
    private HashSet<String> codons;

    /**
     * Construct a new codon filter.
     *
     * @param codon		one or more codons considered acceptable
     */
    public CodonFilter(String... codon) {
        this.codons = new HashSet<String>(codon.length);
        for (String filterCodon : codon) {
            this.codons.add(filterCodon.toUpperCase());
        }
    }

    /**
     * @return TRUE if the codon at the specified position matches the filter
     *
     * @param pos		position (1-based) in a sequence
     * @param sequence	DNA sequence to check
     */
    public boolean matches(int pos, String sequence) {
        String codon = StringUtils.substring(sequence, pos-1, pos+2).toUpperCase();
        return this.codons.contains(codon);
    }

}
