package uk.co.norphos.crystallography.toolkit;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class UnitCell {
	
	private Lattice real, reciprocal;
	private Double volume;
	private RealMatrix metricTensor, reciprocalMetricTensor;

	public UnitCell(Lattice realSpaceLattice) {
		updateCell(realSpaceLattice);
	}
	
	public void updateCell(Lattice realSpaceLattice) {
		this.real = realSpaceLattice;
		volume = null;
	
		metricTensor = determineMetricTensor();
		reciprocalMetricTensor = new LUDecomposition(metricTensor).getSolver().getInverse();
		volume = calculateCellVolume();
	}
	
	private RealMatrix determineMetricTensor() {
		double p00 = real.a()*real.a();
		double p11 = real.b()*real.b();
		double p22 = real.c()*real.c();
		
		double p01 = offAxisCalculator(real.a(), real.b(), real.gaR());
		double p02 = offAxisCalculator(real.a(), real.c(), real.beR());
		double p12 = offAxisCalculator(real.b(), real.c(), real.alR());
		
		return MatrixUtils.createRealMatrix(new double[][]{
			{p00, p01, p02},
			{p01, p11, p12},
			{p02, p12, p22}});
		}
	
	private double offAxisCalculator(double a, double b, double angle) {
		double result = a *b * Math.cos(angle);
		if (Math.abs(result) < 1e-10) return 0.0;
		return result;
	}
	
	private double calculateCellVolume() {
		return Math.sqrt(new LUDecomposition(metricTensor).getDeterminant());
	}
	
	public RealMatrix getMetricTensor() {
		return metricTensor;
	}
	
	public double getCellVolume() {
		return volume;
	}
	
	public double findVectorMagnitude(RealVector vector) {
		return findVectorMagnitude(vector, metricTensor);
	}
	
	public double findVectorMagnitude(RealVector vector, RealMatrix tensor) {
		double product = vector.dotProduct(tensor.operate(vector));
		return Math.sqrt(product);
	}

}
