/**
 *
 */
package org.theseed.genome.contigs;

/**
 * In this class, the sensor value is the three nucleotides at each position before and after the center. The left width
 * must be a multiple of 3 and the right width one less than a multiple of 3.  The triples are all converted to lower case,
 * and anything other than the four real ones is converted to "n".
 *
 * @author Bruce Parrello
 *
 */
public class CodonContigSensorFactory extends ContigSensorFactory {

    @Override
    protected void convertSequence(ContigSensor sensor, String sequence, int pos) {
        boolean suspicion = false;
        int offset = pos - ContigSensorFactory.getLeftWidth() - 1;
        int stride = this.getStride();
        int fullWidth = ContigSensorFactory.getFullWidth();
        String[] buffer = new String[fullWidth / stride];
        for (int i = 0; i < buffer.length; i++) {
            StringBuffer v = new StringBuffer(3);
            for (int j = offset; j <= offset + 2; j++) {
                if (j < 0 || j >= sequence.length()) {
                    v.append('-');
                } else {
                    switch (sequence.charAt(j)) {
                    case 'A' :
                    case 'a' :
                        v.append('a');
                        break;
                    case 'C' :
                    case 'c' :
                        v.append('c');
                        break;
                    case 'G' :
                    case 'g' :
                        v.append('g');
                        break;
                    case 'T' :
                    case 't' :
                    case 'U' :
                    case 'u' :
                        v.append('t');
                        break;
                    default :
                        v.append('n');
                        suspicion = true;
                    }
                }
            }
            buffer[i] = v.toString();
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
