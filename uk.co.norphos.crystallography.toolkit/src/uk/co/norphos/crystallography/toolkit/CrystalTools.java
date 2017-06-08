package uk.co.norphos.crystallography.toolkit;

import org.apache.commons.math3.linear.RealMatrix;

public class CrystalTools {
	
	public static double uijToUiso(RealMatrix uijMatrix) {
		//TODO
		return 0d;
	}
	
	public static double uijToBeq(RealMatrix uijMatrix) {
		return uisoToBeq(uijToUiso(uijMatrix));
	}
	
	public static double beqToUiso(double beq) {
		//TODO
		return 0d;
	}
	
	public static double uisoToBeq(double uiso) {
		//TODO
		return 0;
	}
	
	public static double qToD(double q) {
		return 2 * Math.PI / q;
	}
	
	public static double dToQ(double d) {
		return qToD(d); //The functions are equivalent
	}
	
	/**
	 * 
	 * @param d double in Angstrom
	 * @param lambda double in Angstrom
	 * @return double two theta in degrees
	 */
	public static double dToTwoTheta(double d, double lambda) {
		return Math.toDegrees(2 * Math.asin(lambda / (2 * d)));
	}
	
	/**
	 * 
	 * @param tth double in degrees
	 * @param lambda double in Angstrom
	 * @return double d in Angstrom
	 */
	public static double twoThetaToD(double tth, double lambda) {
		return lambda / (2 * Math.sin(Math.toRadians(tth) / 2));
	}
	
	/**
	 * 
	 * @param q double (Q) in Angstrom<sup>-1</sup>
	 * @param lambda double in Angstrom
	 * @return double two theta in degrees
	 */
	public static double qTotwoTheta(double q, double lambda) {
		return dToTwoTheta(qToD(q), lambda);
	}
	
	/**
	 * 
	 * @param tth double in degrees
	 * @param lambda double in Angstrom
	 * @return double Q in Angstrom<sup>-1</sup>
	 */
	public static double twoThetaToQ(double tth, double lambda) {
		return dToQ(twoThetaToD(tth, lambda));
	}

}
