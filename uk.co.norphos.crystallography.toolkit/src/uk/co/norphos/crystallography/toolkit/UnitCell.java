package uk.co.norphos.crystallography.toolkit;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class UnitCell {
	
	private Lattice real, reciprocal;
	private Double volume, reciprocalVolume;
	private RealMatrix metricTensor, reciprocalMetricTensor;
	private LUDecomposition metricTensorLUDecomp, reciprocalMetricTensorLUDecomp;

	public UnitCell(Lattice realSpaceLattice) {
		updateCell(realSpaceLattice);
	}
	
	public void updateCell(Lattice realSpaceLattice) {
		this.real = realSpaceLattice;
		volume = null;
	
		metricTensor = determineMetricTensor();
		metricTensorLUDecomp = new LUDecomposition(metricTensor);
		volume = Math.sqrt(metricTensorLUDecomp.getDeterminant());
		reciprocalMetricTensor = metricTensorLUDecomp.getSolver().getInverse();
		reciprocalMetricTensorLUDecomp = new LUDecomposition(reciprocalMetricTensor);
		reciprocalVolume = Math.sqrt(reciprocalMetricTensorLUDecomp.getDeterminant());
		reciprocal = determineReciprocalLattice();
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
	
	private Lattice determineReciprocalLattice() {
		double rA = Math.sqrt(reciprocalMetricTensor.getEntry(0, 0));
		double rB = Math.sqrt(reciprocalMetricTensor.getEntry(1, 1));
		double rC = Math.sqrt(reciprocalMetricTensor.getEntry(2, 2));
		double rAl = Math.toDegrees(Math.acos(reciprocalMetricTensor.getEntry(1, 2) / (rB * rC)));
		double rBe = Math.toDegrees(Math.acos(reciprocalMetricTensor.getEntry(0, 2) / (rA * rC)));
		double rGa = Math.toDegrees(Math.acos(reciprocalMetricTensor.getEntry(0, 1) / (rA * rB)));
		
		return new Lattice.LatticeBuilder(rA).setB(rB).setC(rC).setAl(rAl).setBe(rBe).setGa(rGa).build();
	}
	
	private double offAxisCalculator(double a, double b, double angle) {
		double result = a *b * Math.cos(angle);
		if (Math.abs(result) < 1e-10) return 0.0;
		return result;
	}
	
	public Lattice getCell() {
		return real;
	}
	
	public RealMatrix getMetricTensor() {
		return metricTensor;
	}
	
	public Lattice getReciprocalLattice() {
		return reciprocal;
	}
	
	public RealMatrix getReciprocalMetricTensor() {
		return reciprocalMetricTensor;
	}
	
	public double getCellVolume() {
		return volume;
	}
	
	public double getReciprocalCellVolume() {
		return reciprocalVolume;
	}
	
	public double findVectorMagnitude(RealVector vector) {
		return findVectorMagnitude(vector, metricTensor);
	}
	
	public double findPlaneDSpacing(RealVector hklVector) {
		return 1/findVectorMagnitude(hklVector, reciprocalMetricTensor);
	}
	
	private double findVectorMagnitude(RealVector vector, RealMatrix tensor) {
		double product = vector.dotProduct(tensor.operate(vector));
		return Math.sqrt(product);
	}

}
