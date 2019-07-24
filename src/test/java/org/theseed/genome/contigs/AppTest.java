package org.theseed.genome.contigs;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.List;
import org.theseed.genomes.Contig;
import org.theseed.locations.Frame;
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
        // Test plus-only mode.
        ContigSensor.setPlusOnly(true);
        assertThat("Plus-only did not adjust width.", ContigSensor.getSensorWidth(), equalTo(5));
        ByteArrayOutputStream testStream = new ByteArrayOutputStream();
        PrintStream testPrinter = new PrintStream(testStream);
        ContigSensor.sensor_headers(testPrinter);
        assertThat("Wrong half-header", testStream.toString(), equalTo("\tpos.00\tpos.01\tpos.02\tpos.03\tpos.04"));
        sensor1 = new ContigSensor("ABC", 2, "ACGTTAGGTT");
        assertThat("Wrong half-sensors", sensor1.getSensorList(),
                contains(-0.3, 0.3, 0.6, 0.6, -0.6));
        assertThat("Wrong string", sensor1.toString(), equalTo("-0.3\t0.3\t0.6\t0.6\t-0.6"));
    }

    /*
     * Test location sensors.
     */
    public void testLocationClasses() {
    	LocationClass lsensor = LocationClass.scheme(LocationClass.Type.CODING, true);
    	assertThat("Incorrect normalize for M2, negative.", lsensor.normalize(Frame.M2), equalTo(Frame.M2));
    	assertThat("Incorrect normalize for M1, negative.", lsensor.normalize(Frame.M1), equalTo(Frame.M1));
    	assertThat("Incorrect normalize for M0, negative.", lsensor.normalize(Frame.M0), equalTo(Frame.M0));
    	assertThat("Incorrect normalize for F0, negative.", lsensor.normalize(Frame.F0), equalTo(Frame.F0));
    	assertThat("Incorrect normalize for P2, negative.", lsensor.normalize(Frame.P2), equalTo(Frame.P2));
    	assertThat("Incorrect normalize for P1, negative.", lsensor.normalize(Frame.P1), equalTo(Frame.P1));
    	assertThat("Incorrect normalize for P0, negative.", lsensor.normalize(Frame.P0), equalTo(Frame.P0));
    	assertThat("Incorrect normalize for XX, negative.", lsensor.normalize(Frame.XX), equalTo(Frame.XX));
    	assertThat("Incorrect coding result for M2, negative.", lsensor.classOf(Frame.P0, Frame.M2), equalTo("coding"));
    	assertThat("Incorrect coding result for M1, negative.", lsensor.classOf(Frame.P0, Frame.M1), equalTo("coding"));
    	assertThat("Incorrect coding result for M0, negative.", lsensor.classOf(Frame.P0, Frame.M0), equalTo("coding"));
    	assertThat("Incorrect coding result for F0, negative.", lsensor.classOf(Frame.P0, Frame.F0), equalTo("space"));
    	assertThat("Incorrect coding result for P2, negative.", lsensor.classOf(Frame.P0, Frame.P2), equalTo("coding"));
    	assertThat("Incorrect coding result for P1, negative.", lsensor.classOf(Frame.P0, Frame.P1), equalTo("coding"));
    	assertThat("Incorrect coding result for P0, negative.", lsensor.classOf(Frame.P0, Frame.P0), equalTo("coding"));
    	assertNull("Incorrect coding result for XX, negative.", lsensor.classOf(Frame.P0, Frame.XX));
    	lsensor = LocationClass.scheme(LocationClass.Type.CODING, false);
    	assertThat("Incorrect coding result for M2, positive.", lsensor.classOf(Frame.P0, Frame.M2), equalTo("space"));
    	assertThat("Incorrect coding result for M1, positive.", lsensor.classOf(Frame.P0, Frame.M1), equalTo("space"));
    	assertThat("Incorrect coding result for M0, positive.", lsensor.classOf(Frame.P0, Frame.M0), equalTo("space"));
    	assertThat("Incorrect coding result for F0, positive.", lsensor.classOf(Frame.P0, Frame.F0), equalTo("space"));
    	assertThat("Incorrect coding result for P2, positive.", lsensor.classOf(Frame.P0, Frame.P2), equalTo("coding"));
    	assertThat("Incorrect coding result for P1, positive.", lsensor.classOf(Frame.P0, Frame.P1), equalTo("coding"));
    	assertThat("Incorrect coding result for P0, positive.", lsensor.classOf(Frame.P0, Frame.P0), equalTo("coding"));
    	assertThat("Incorrect normalize for M2, positive.", lsensor.normalize(Frame.M2), equalTo(Frame.F0));
    	assertThat("Incorrect normalize for M1, positive.", lsensor.normalize(Frame.M1), equalTo(Frame.F0));
    	assertThat("Incorrect normalize for M0, positive.", lsensor.normalize(Frame.M0), equalTo(Frame.F0));
    	assertThat("Incorrect normalize for F0, positive.", lsensor.normalize(Frame.F0), equalTo(Frame.F0));
    	assertThat("Incorrect normalize for P2, positive.", lsensor.normalize(Frame.P2), equalTo(Frame.P2));
    	assertThat("Incorrect normalize for P1, positive.", lsensor.normalize(Frame.P1), equalTo(Frame.P1));
    	assertThat("Incorrect normalize for P0, positive.", lsensor.normalize(Frame.P0), equalTo(Frame.P0));
    	assertThat("Incorrect normalize for XX, positive.", lsensor.normalize(Frame.XX), equalTo(Frame.XX));
    	assertNull("Incorrect coding result for XX, positive.", lsensor.classOf(Frame.P0, Frame.XX));
    	lsensor = LocationClass.scheme(LocationClass.Type.PHASE, true);
    	assertThat("Incorrect coding result for M2, negative.", lsensor.classOf(Frame.P0, Frame.M2), equalTo("-3"));
    	assertThat("Incorrect coding result for M1, negative.", lsensor.classOf(Frame.P0, Frame.M1), equalTo("-2"));
    	assertThat("Incorrect coding result for M0, negative.", lsensor.classOf(Frame.P0, Frame.M0), equalTo("-1"));
    	assertThat("Incorrect coding result for F0, negative.", lsensor.classOf(Frame.P0, Frame.F0), equalTo("0"));
    	assertThat("Incorrect coding result for P2, negative.", lsensor.classOf(Frame.P0, Frame.P2), equalTo("+3"));
    	assertThat("Incorrect coding result for P1, negative.", lsensor.classOf(Frame.P0, Frame.P1), equalTo("+2"));
    	assertThat("Incorrect coding result for P0, negative.", lsensor.classOf(Frame.P0, Frame.P0), equalTo("+1"));
    	assertNull("Incorrect coding result for XX, negative.", lsensor.classOf(Frame.P0, Frame.XX));
    	lsensor = LocationClass.scheme(LocationClass.Type.PHASE, false);
    	assertThat("Incorrect coding result for M2, positive.", lsensor.classOf(Frame.P0, Frame.M2), equalTo("0"));
    	assertThat("Incorrect coding result for M1, positive.", lsensor.classOf(Frame.P0, Frame.M1), equalTo("0"));
    	assertThat("Incorrect coding result for M0, positive.", lsensor.classOf(Frame.P0, Frame.M0), equalTo("0"));
    	assertThat("Incorrect coding result for F0, positive.", lsensor.classOf(Frame.P0, Frame.F0), equalTo("0"));
    	assertThat("Incorrect coding result for P2, positive.", lsensor.classOf(Frame.P0, Frame.P2), equalTo("+3"));
    	assertThat("Incorrect coding result for P1, positive.", lsensor.classOf(Frame.P0, Frame.P1), equalTo("+2"));
    	assertThat("Incorrect coding result for P0, positive.", lsensor.classOf(Frame.P0, Frame.P0), equalTo("+1"));
       	assertNull("Incorrect coding result for XX, positive.", lsensor.classOf(Frame.P0, Frame.XX));
           	try {
    		lsensor = LocationClass.scheme(LocationClass.Type.EDGE, true);
    		fail("Invalid combination passed.");
    	} catch (IllegalArgumentException e) { }
    	lsensor = LocationClass.scheme(LocationClass.Type.EDGE, false);
    	assertThat("Incorrect coding result for M2, positive.", lsensor.classOf(Frame.F0, Frame.M2), equalTo("other"));
    	assertThat("Incorrect coding result for M1, positive.", lsensor.classOf(Frame.M2, Frame.M1), equalTo("other"));
    	assertThat("Incorrect coding result for M0, positive.", lsensor.classOf(Frame.M1, Frame.M0), equalTo("other"));
    	assertThat("Incorrect coding result for M0/M2, positive.", lsensor.classOf(Frame.M0, Frame.M2), equalTo("other"));
    	assertThat("Incorrect coding result for M0/F0, positive.", lsensor.classOf(Frame.M0, Frame.F0), equalTo("other"));
    	assertThat("Incorrect coding result for P0, positive.", lsensor.classOf(Frame.F0, Frame.P0), equalTo("start"));
    	assertThat("Incorrect coding result for P1, positive.", lsensor.classOf(Frame.P0, Frame.P1), equalTo("other"));
    	assertThat("Incorrect coding result for P2, positive.", lsensor.classOf(Frame.P1, Frame.P2), equalTo("other"));
       	assertThat("Incorrect coding result for P2/P0, positive.", lsensor.classOf(Frame.P2, Frame.P0), equalTo("other"));
     	assertThat("Incorrect coding result for P2/F0, positive.", lsensor.classOf(Frame.P2, Frame.F0), equalTo("stop"));
       	assertNull("Incorrect coding result for XX, positive.", lsensor.classOf(Frame.P0, Frame.XX));

    }
}
