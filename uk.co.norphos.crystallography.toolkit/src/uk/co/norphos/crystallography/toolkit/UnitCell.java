package uk.co.norphos.crystallography.toolkit;

import org.eclipse.january.dataset.Dataset;
import org.eclipse.january.dataset.DatasetFactory;
import org.eclipse.january.dataset.LinearAlgebra;

public class UnitCell {
	
	private Lattice real, reciprocal;
	private Double volume;
	private Dataset metricTensor, reciprocalMetricTensor;

	public UnitCell(Lattice realSpaceLattice) {
		updateCell(realSpaceLattice);
	}
	
	public void updateCell(Lattice realSpaceLattice) {
		this.real = realSpaceLattice;
		volume = null;
	
		metricTensor = determineMetricTensors();
		reciprocalMetricTensor = LinearAlgebra.calcInverse(metricTensor);
		volume = calculateCellVolume();
	}
	
	private Dataset determineMetricTensors() {
		Dataset metricTensor = DatasetFactory.zeros(new int[]{3, 3});
		
		double p00 = real.a()*real.a();
		double p11 = real.b()*real.b();
		double p22 = real.c()*real.c();
		
		double p01 = offAxisCalculator(real.a(), real.b(), real.gaR());
		double p02 = offAxisCalculator(real.a(), real.c(), real.beR());
		double p12 = offAxisCalculator(real.b(), real.c(), real.alR());
		
		metricTensor.set(p00, 0, 0);
		metricTensor.set(p01, 0, 1);
		metricTensor.set(p02, 0, 2);
		metricTensor.set(p01, 1, 0);
		metricTensor.set(p11, 1, 1);
		metricTensor.set(p12, 1, 2);
		metricTensor.set(p02, 2, 0);
		metricTensor.set(p12, 2, 1);
		metricTensor.set(p22, 2, 2);
		
		return metricTensor;
	}
	
	private double offAxisCalculator(double a, double b, double angle) {
		double result = a *b * Math.cos(angle);
		if (Math.abs(result) < 1e-10) return 0.0;
		return result;
	}
	
	private double calculateCellVolume() {
		return Math.sqrt(LinearAlgebra.calcDeterminant(metricTensor));
	}
	
	public Dataset getMetricTensor() {
		return metricTensor;
	}
	
	public double getCellVolume() {
		return volume;
	}

}
