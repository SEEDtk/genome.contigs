/**
 *
 */
package org.theseed.genome.contigs;

import java.util.HashMap;

/**
 * In this class, the sensor value is the amino acid at each position before and after the center. The left width
 * must be a multiple of 3 and the right width one less than a multiple of 3.  The amino acids are all converted
 * to upper case.  Anything off the edge of the contig converts to "-" and anything with an ambiguity character
 * converts to "X".  Currently, only genetic code 11 is supported.
 *
 * @author Bruce Parrello
 *
 */
public class AminoAcidContigSensorFactory extends ContigSensorFactory {

    /** genetic code translation table */
    @SuppressWarnings("serial")
    private static final HashMap<String, String> GENETIC_CODE_11 = new HashMap<String, String>() {{
        put("AAA","K"); put("AAC","N"); put("AAG","K"); put("AAT","N"); put("ACA","T");
        put("ACC","T"); put("ACG","T"); put("ACT","T"); put("AGA","R"); put("AGC","S");
        put("AGG","R"); put("AGT","S"); put("ATA","I"); put("ATC","I"); put("ATG","M");
        put("ATT","I"); put("CAA","Q"); put("CAC","H"); put("CAG","Q"); put("CAT","H");
        put("CCA","P"); put("CCC","P"); put("CCG","P"); put("CCT","P"); put("CGA","R");
        put("CGC","R"); put("CGG","R"); put("CGT","R"); put("CTA","L"); put("CTC","L");
        put("CTG","L"); put("CTT","L"); put("GAA","E"); put("GAC","D"); put("GAG","E");
        put("GAT","D"); put("GCA","A"); put("GCC","A"); put("GCG","A"); put("GCT","A");
        put("GGA","G"); put("GGC","G"); put("GGG","G"); put("GGT","G"); put("GTA","V");
        put("GTC","V"); put("GTG","V"); put("GTT","V"); put("TAA","*"); put("TAC","Y");
        put("TAG","*"); put("TAT","Y"); put("TCA","S"); put("TCC","S"); put("TCG","S");
        put("TCT","S"); put("TGA","*"); put("TGC","C"); put("TGG","W"); put("TGT","C");
        put("TTA","L"); put("TTC","F"); put("TTG","L"); put("TTT","F");
    }};


    @Override
    protected void convertSequence(ContigSensor sensor, String sequence, int pos) {
        boolean suspicion = false;
        int offset = pos - ContigSensorFactory.getLeftWidth() - 1;
        int stride = this.getStride();
        int fullWidth = ContigSensorFactory.getFullWidth();
        String[] buffer = new String[fullWidth / stride];
        for (int i = 0; i < buffer.length; i++) {
            String aa = "-";
            int endOffset = offset + 3;
            if (offset >= 0 && endOffset <= sequence.length()) {
                String codon = sequence.substring(offset, endOffset).toUpperCase();
                aa = GENETIC_CODE_11.get(codon);
                if (aa == null) {
                    aa = "X";
                    suspicion = true;
                }
            }
            buffer[i] = aa;
            offset += stride;
        }
        sensor.storeSensors(buffer, suspicion);
    }

    /**
     * @return the stride between positions to be sensed (3 for this type)
     */
    protected int getStride() {
        return 3;
    }

}
