package uk.co.norphos.crystallography.toolkit;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Before;
import org.junit.Test;

import uk.co.norphos.crystallography.api.CrystalSystem;
import uk.co.norphos.crystallography.api.Lattice;

public class CrystallographyFactoryTest {
	
	private Lattice cubicLatt, tricLatt;
	
	@Before
	public void setUp() {
		/*
		 * Structures from which data taken are located in the test-data directory
		 * Volumes calculated using VESTA
		 *///cubicVol = 160.11892937380182
		cubicLatt =  new Lattice(5.43018, 5.43018, 5.43018, 90, 90, 90, 160.118936, CrystalSystem.CUBIC);
		tricLatt = new Lattice(7.19196, 8.12720, 8.12771, 82.4809, 69.2610, 69.2584, 415.482298, CrystalSystem.TRICLINIC);
	}

	@Test
	public void testCreateLatticeFromTensor() {
		RealMatrix cubicTensor = makeMetricTensor(cubicLatt);
		Lattice madeLattice = CrystallographyFactory.createLattice(cubicTensor);
		assertEquals(cubicLatt, madeLattice);
		
		cubicTensor = makeMetricTensor(tricLatt);
		madeLattice = CrystallographyFactory.createLattice(cubicTensor);
		assertEquals(tricLatt, madeLattice);
	}
	
	private RealMatrix makeMetricTensor(Lattice lattice) {
		double[] lengths = lattice.getLengths();
		double[] angles = lattice.getAnglesRadians();
		
		double[][] tensor = new double[3][3];
		int i, j, k;
		//Set diagonal values
		for (i = 0; i < 3; i++) {
			tensor[i][i] = Math.pow(lengths[i], 2);
		}
		//Set off diagonal values
		for (i = 0; i < 3; i++) {
			j = (i + 1) % 3;
			k = (j + 1) % 3;
			double val = lengths[i] * lengths[j] * Math.cos(angles[k]);
			if (Math.abs(val) < 1e-10) val = 0;
			tensor[i][j] = val;
			tensor[j][i] = val;
		}
		return MatrixUtils.createRealMatrix(tensor);
	}

}
