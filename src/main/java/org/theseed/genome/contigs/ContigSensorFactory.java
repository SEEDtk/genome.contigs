/**
 *
 */
package org.theseed.genome.contigs;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genomes.Contig;
import org.theseed.sequence.Sequence;

/**
 * @author Bruce Parrello
 *
 */
public abstract class ContigSensorFactory {

    /**
     * types of sensors currently supported
     */
    public static enum Type {
        DIRECT, ONE_HOT, CODON
    }

    /** global sensor width, to either side of the target position */
    protected static int halfWidth = 14;

    /**
     * Construct a blank, empty sensor factory.
     */
    public ContigSensorFactory() {
    }

    /**
     * @return the number of sensors on each side of the target position
     */
    public static int getHalfWidth() {
        return halfWidth;
    }

    /**
     * Specify a new global half-width.  This does not affect sensors
     * already constructed.
     *
     * @param new sensor width
     */
    public static void setHalfWidth(int newWidth) {
        ContigSensorFactory.halfWidth = newWidth;
    }

    /**
     * @return a list of all the non-suspicious sensors for the specified sequence
     *
     * @param inSeq		DNA sequence to process
     */
    public List<ContigSensor> processContig(Sequence inSeq) {
        return processContig(inSeq.getLabel(), inSeq.getSequence(), 1, inSeq.length());
    }

    /**
     * @return a list of all the non-suspicious sensors for the specified portion of a DNA sequence
     *
     * @param contigId	contig ID of the sequence
     * @param sequence	DNA of the sequence
     * @param start		starting position
     * @param len		number of positions to process
     */
    private List<ContigSensor> processContig(String contigId, String sequence, int start, int len) {
        ArrayList<ContigSensor> retVal = new ArrayList<ContigSensor>(sequence.length());
        int end = start + len - 1;
        if (end > sequence.length()) end = sequence.length();
        for (int i = start; i <= end; i++) {
            ContigSensor snapshot = this.create(contigId, i, sequence);
            if (! snapshot.suspicious) retVal.add(snapshot);
        }
        return retVal;
    }

    /**
     * Create a contig sensor of the appropriate type at the specified position in a DNA sequence
     *
     * @param id		ID of the DNA sequence
     * @param pos		position in the sequence for the sensor
     * @param sequence	DNA sequence from which the sensor is derived
     */
    public ContigSensor create(String id, int pos, String sequence) {
        ContigSensor retVal = new ContigSensor(id, pos);
        this.convertSequence(retVal, sequence, pos);
        return retVal;
    }

    /**
     * Fill the specified sensor with the DNA information at the current position.
     *
     * @param sensor	target contig sensor
     * @param sequence	source sequence being sensed
     * @param pos		position of the sensor
     */
    protected abstract void convertSequence(ContigSensor sensor, String sequence, int pos);

    /**
     * @return a list of the non-suspicious sensors for the specified contig region
     *
     * @param contig	source contig
     * @param start		position for the first sensor
     * @param len		number of positions to process
     */
    public List<ContigSensor> processContig(Contig contig, int start, int len) {
        return processContig(contig.getId(), contig.getSequence(), start, len);
    }

    /**
     * @return a sensor factory of the specified type
     *
     * @param type	type of contig sensors to create
     */
    public static ContigSensorFactory create(Type type) {
        ContigSensorFactory retVal = null;
        switch (type) {
        case DIRECT :
            retVal = new DirectContigSensorFactory();
            break;
        case ONE_HOT :
            retVal = new OneHotContigSensorFactory();
            break;
        case CODON :
            retVal = new CodonContigSensorFactory();
            break;
        default :
            throw new IllegalArgumentException("Unknown contig factory type " + type + ".");
        }
        return retVal;
    }

    /**
     * @return the sensor column headers for this sensor type
     */
    public String sensor_headers() {
        int halfWidth = ContigSensorFactory.getHalfWidth();
        int fullWidth = halfWidth * 2 + 1;
        ArrayList<String> headers = new ArrayList<String>(fullWidth);
        for (int i = -halfWidth; i <= halfWidth; i++) {
            headers.add("pos." + i);
        }
        return StringUtils.join(headers, '\t');
    }

}
