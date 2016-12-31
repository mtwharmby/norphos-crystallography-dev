package uk.co.norphos.crystallography.toolkit;

import static org.junit.Assert.assertEquals;

import org.eclipse.january.asserts.TestUtils;
import org.eclipse.january.dataset.Dataset;
import org.eclipse.january.dataset.DatasetFactory;
import org.junit.Before;
import org.junit.Test;

public class UnitCellTest {
	
	private UnitCell uc;
	
	@Before
	public void setUp() {
		
	}
	
	@Test
	public void testMetricTensor() {
		//Cubic a = 3
		uc = new UnitCell(new Lattice.LatticeBuilder(3).build());
		Dataset fakeTensor = DatasetFactory.createFromObject(new double[]{9,0,0, 0,9,0, 0,0,9}, 3,3);
		TestUtils.assertDatasetEquals(fakeTensor, uc.getMetricTensor(), 0, 1e-10);
		
		//TODO Rhombohedral a = 3; al = ?
		//TODO Hexagonal a = 5; c = 2; ga = 120
		
		//Tetragonal a = b = 2; c = 5
		uc.updateCell(new Lattice.LatticeBuilder(2).setC(5).build());
		fakeTensor = DatasetFactory.createFromObject(new double[]{4,0,0, 0,4,0, 0,0,25}, 3,3);
		TestUtils.assertDatasetEquals(fakeTensor, uc.getMetricTensor(), 0, 1e-10);
		
		//Orthorhombic a = 5; b = 3; c = 2
		uc.updateCell(new Lattice.LatticeBuilder(5).setB(3).setC(2).build());
		fakeTensor = DatasetFactory.createFromObject(new double[]{25,0,0, 0,9,0, 0,0,4}, 3,3);
		TestUtils.assertDatasetEquals(fakeTensor, uc.getMetricTensor(), 0, 1e-10);
		
		//Monoclinic a = 2; b = 3; c = 5; be = 30deg
		uc.updateCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).setBe(30).build());
		fakeTensor = DatasetFactory.createFromObject(new double[]
					{4,0,5*Math.sqrt(3), 
					 0,9,0, 
					 5*Math.sqrt(3),0,25}, 
				3,3);
		TestUtils.assertDatasetEquals(fakeTensor, uc.getMetricTensor(), 0, 1e-10);
		
		//Triclinic #a = 3; b = 5; c = 2; al = 30deg; be = 45 deg; ga = 60 deg
		uc.updateCell(new Lattice.LatticeBuilder(3).setB(5).setC(2).setAl(30).setBe(45).setGa(60).build());
		fakeTensor = DatasetFactory.createFromObject(new double[]{
					9,7.5,6/Math.sqrt(2),
					7.5, 25, 5*Math.sqrt(3),
					6/Math.sqrt(2), 5*Math.sqrt(3), 4},
				3,3);
		TestUtils.assertDatasetEquals(fakeTensor, uc.getMetricTensor(), 0, 1e-10);
	}
	
	@Test
	public void testCellVolume() {
		uc = new UnitCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).build());
		assertEquals(30.0, uc.getCellVolume(), 0.05);
	}

}
