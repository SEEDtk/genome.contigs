package org.theseed.genome.contigs;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.util.List;
import org.theseed.genomes.Contig;
import org.theseed.sequence.Sequence;



/**
 * Unit test for simple App.
 */
public class AppTest
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }


    /**
     * Test the contig sensors.
     */
    public void testContigSensors()
    {
        String contigID = "3000.contig.1";
        Sequence frec = new Sequence(contigID, "", "AACGTCCTGAAGTC");
        ContigSensor.setHalfWidth(4);
        assertThat("Wrong computed width.", ContigSensor.getSensorWidth(), equalTo(9));
        List<ContigSensor> sensors = ContigSensor.processContig(frec);
        ContigSensor sensor1 = sensors.get(0);
        ContigSensor sensor11 = sensors.get(11);
        assertThat("Wrong contig ID for first sensor.", sensor1.getContigId(), equalTo(contigID));
        assertThat("Wrong position for first sensor.", sensor1.getPosition(), equalTo(1));
        assertThat("Wrong meta string for first sensor.", sensor1.getMeta(), equalTo(contigID + ";1"));
        assertThat("Wrong sensors for first", sensor1.getSensorList(),
                contains(0.0, 0.0, 0.0, 0.0, -0.6, -0.6, -0.3, 0.3, 0.6));
        assertThat("Wrong sensors for 11th", sensor11.getSensorList(),
                contains(0.6, 0.3, -0.6, -0.6, 0.3, 0.6, -0.3, 0.0, 0.0));
        assertThat("Wrong contig ID for 11th", sensor11.getContigId(), equalTo(contigID));
        assertThat("Wrong position for 11th", sensor11.getPosition(), equalTo(12));
        for (ContigSensor sensor : sensors) {
            assertThat("Wrong contig ID in sensor", sensor.getContigId(), equalTo(contigID));
        }
        Contig contig = new Contig(contigID, "AACGTNCGGGGAAAT", 11);
        sensors = ContigSensor.processContig(contig, 2, 20);
        sensor11 = sensors.get(0);
        assertThat("Wrong position for post-N sensor.", sensor11.getPosition(), equalTo(11));
        sensor1 = new ContigSensor(contigID, 7, contig.getSequence());
        assertTrue("Sensor not suspicious.", sensor1.isSuspicious());
        double[] sensorArray = sensor1.getSensors();
        List<Double> sensorList = sensor1.getSensorList();
        for (int i = 0; i < ContigSensor.getSensorWidth(); i++) {
            assertThat("Sensor mismatch at " + i, sensorArray[i], equalTo(sensorList.get(i)));
        }
    }
}
