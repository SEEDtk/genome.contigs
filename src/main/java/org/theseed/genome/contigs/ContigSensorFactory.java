/**
 *
 */
package org.theseed.genome.contigs;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Contig;
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
        DIRECT, CHANNEL, CODON, AMINOACID
    }

    /** global sensor width, to either side of the target position */
    protected static int leftWidth = 21;
    protected static int rightWidth = 45;

    /**
     * Construct a blank, empty sensor factory.
     */
    public ContigSensorFactory() {
    }

    /**
     * @return the number of sensors on the left side of the target position
     */
    public static int getLeftWidth() {
        return leftWidth;
    }

    /**
     * @return the number of sensors on the right side of the target position
     */
    public static int getRightWidth() {
        return rightWidth;
    }

    /**
     * @return the full width of a sensor
     */
    public static int getFullWidth() {
        return leftWidth + rightWidth + 1;
    }

    /**
     * Specify a new global left-width.  This does not affect sensors
     * already constructed.
     *
     * @param new sensor width
     */
    public static void setLeftWidth(int newWidth) {
        ContigSensorFactory.leftWidth = newWidth;
    }

    /**
     * Specify a new global right-width.  This does not affect sensors
     * already constructed.
     *
     * @param new sensor width
     */
    public static void setRightWidth(int newWidth) {
        ContigSensorFactory.rightWidth = newWidth;
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
     * @return the stride between positions inside the sensor (normnally 1, sometimes 3)
     */
    protected int getStride() {
        return 1;
    }

    /**
     * Create a contig sensor of the appropriate type at the specified position in a DNA sequence
     *
     * @param id		ID of the DNA sequence
     * @param pos		position in the sequence for the sensor
     * @param sequence	DNA sequence from which the sensor is derived
     */
    public ContigSensor create(String id, int pos, String sequence) {
        String codon = CodonFilter.getCodon(pos, sequence);
        ContigSensor retVal = new ContigSensor(id, pos, codon);
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
        case CHANNEL :
            retVal = new ChannelContigSensorFactory();
            break;
        case CODON :
            retVal = new CodonContigSensorFactory();
            break;
        case AMINOACID :
            retVal = new AminoAcidContigSensorFactory();
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
        ArrayList<String> headers = new ArrayList<String>(getFullWidth());
        int stride = this.getStride();
        for (int i = -getLeftWidth(); i <= getRightWidth(); i += stride) {
            headers.add("pos." + i);
        }
        return StringUtils.join(headers, '\t');
    }

}
