/**
 *
 */
package org.theseed.genome.contigs;

import org.theseed.locations.Frame;

/**
 * This is the base class for location classifiers.  It computes the class of a location
 * (usually related to its coding function) based on the current and previous frame
 * codes.  The constructor specifies whether or not minus strand proteins are considered
 * as coding regions.
 *
 * @author Bruce Parrello
 *
 */
public abstract class LocationClass {

	// FIELDS
	/** If TRUE, then minus-strand proteins are considered coding regions */
	boolean	negative;

	/**
	 * Construct a blank location class handler.
	 *
	 * @param negativeFlag	TRUE if minus-strand proteins are considered coding regions
	 */
	public LocationClass(boolean negativeFlag) {
		this.negative = negativeFlag;
	}

	/**
	 * Return the location class.  This method prepares the inputs and asks the
	 * subclass for the result.  A return of NULL means the current location
	 * is invalid.
	 *
	 * @param prevFrame		the frame of the previous location
	 * @param thisFrame		the frame of the current location
	 */
	public String classOf(Frame prevFrame, Frame thisFrame) {
		String retVal = null;
		if (thisFrame != Frame.XX)
			retVal = computeClass(normalize(prevFrame), normalize(thisFrame));
		return retVal;
	}

	/**
	 * @return the class for a specified location
	 *
	 * @param prevFrame	normalized frame of previous location
	 * @param thisFrame	normalized frame of this location
	 */
	protected abstract String computeClass(Frame prevFrame, Frame thisFrame);

	/**
	 * Convert a frame according to the policy on the minus strand.  This
	 * will either return the original frame or Frame.F0.
	 *
	 * @param frm	frame in question
	 */
	public Frame normalize(Frame frm) {
		Frame retVal = frm;
		if (! this.negative && frm.negative()) {
			retVal = Frame.F0;
		}
		return retVal;
	}

	/** enum used to choose a location class scheme */
	public static enum Type {
		/** identify "coding" or "space" (non-coding) locations */
		CODING,
		/** identify "start" (start of coding), "stop" (past end of coding), or "other" locations */
		EDGE,
		/** identify actual coding frame */
		PHASE
	}

	public static LocationClass scheme(Type type, boolean negativeFlag) {
		LocationClass retVal;
		switch (type) {
		case CODING :
			retVal = new LocationClass.Coding(negativeFlag);
			break;
		case EDGE :
			retVal = new LocationClass.Edge(negativeFlag);
			break;
		case PHASE :
			retVal = new LocationClass.Phase(negativeFlag);
			break;
		default :
			retVal = null;
		}
		return retVal;
	}

	// SUBCLASSES

	/**
	 * Classification is the actual coding frame string (+1, +2, +3, 0)
	 */
	public static class Phase extends LocationClass {

		public Phase(boolean negativeFlag) {
			super(negativeFlag);
		}

		@Override
		protected String computeClass(Frame prevFrame, Frame thisFrame) {
			return thisFrame.toString();
		}

	}

	/**
	 * Classification is "start", "stop", or "other".
	 */
	public static class Edge extends LocationClass {

		public Edge(boolean negativeFlag) {
			super(negativeFlag);
			if (this.negative)
				throw new IllegalArgumentException("Cannot use edge mode with the negative strand turned on.");
		}

		@Override
		protected String computeClass(Frame prevFrame, Frame thisFrame) {
			String retVal = "other";
			if (thisFrame == Frame.P0 && prevFrame == Frame.F0)
				retVal = "start";
			else if (prevFrame == Frame.P2 && thisFrame == Frame.F0)
				retVal = "stop";
			return retVal;
		}

	}

	/**
	 * Classification is "coding" or "space".
	 */
	public static class Coding extends LocationClass {

		public Coding(boolean negativeFlag) {
			super(negativeFlag);
		}

		@Override
		protected String computeClass(Frame prevFrame, Frame thisFrame) {
			return (thisFrame != Frame.F0 ? "coding" : "space");
		}

	}


}
