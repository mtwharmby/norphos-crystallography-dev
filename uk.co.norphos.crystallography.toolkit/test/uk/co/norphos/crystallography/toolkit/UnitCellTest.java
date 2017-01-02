package uk.co.norphos.crystallography.toolkit;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.TestUtils;
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
		RealMatrix fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{9,0,0},{0,9,0},{0,0,9}});
		TestUtils.assertEquals("Cubic tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
		
		//Rhombohedral a = 3; al = 60
		uc.updateCell(new Lattice.LatticeBuilder(3).setAl(60).build());
		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{9,4.5,4.5},{4.5,9,4.5},{4.5,4.5,9}});
		TestUtils.assertEquals("Rhombohedral tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);

		//Hexagonal a = 5; c = 2; ga = 120
		uc.updateCell(new Lattice.LatticeBuilder(5).setC(2).setGa(120).build());
		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{25,-12.5,0},{-12.5,25,0},{0,0,4}});
		TestUtils.assertEquals("Hexagonal tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
		
		//Tetragonal a = b = 2; c = 5
		uc.updateCell(new Lattice.LatticeBuilder(2).setC(5).build());
		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{4,0,0},{0,4,0},{0,0,25}});
		TestUtils.assertEquals("Tetragonal tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
		
		//Orthorhombic a = 5; b = 3; c = 2
		uc.updateCell(new Lattice.LatticeBuilder(5).setB(3).setC(2).build());
		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{25,0,0},{0,9,0},{0,0,4}});
		TestUtils.assertEquals("Orthorhombic  tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
		
		//Monoclinic a = 2; b = 3; c = 5; be = 30deg
		uc.updateCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).setBe(30).build());
		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{4,0,5*Math.sqrt(3)},{0,9,0},{5*Math.sqrt(3),0,25}});
		TestUtils.assertEquals("Monoclinic tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
		
		//Triclinic #a = 3; b = 5; c = 2; al = 30deg; be = 45 deg; ga = 60 deg
		uc.updateCell(new Lattice.LatticeBuilder(3).setB(5).setC(2).setAl(30).setBe(45).setGa(60).build());
		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{9,7.5,6/Math.sqrt(2)},{7.5, 25, 5*Math.sqrt(3)},{6/Math.sqrt(2), 5*Math.sqrt(3), 4}});
		TestUtils.assertEquals("Triclinic tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
	}
	
	@Test
	public void testCellVolume() {
		uc = new UnitCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).build());
		assertEquals(30.0, uc.getCellVolume(), 0.05);
	}
	
	@Test
	public void testVectorMagnitude() {
		uc = new UnitCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).build());
		assertEquals("Orthorhombic [100] vector incorrectly calculated", 2.0, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,0,0})), 1e-6);
		assertEquals("Orthorhombic [010] vector incorrectly calculated", 3.0, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{0,1,0})), 1e-6);
		assertEquals("Orthorhombic [110] vector incorrectly calculated", Math.sqrt(4+9), uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,0})), 1e-6);
		assertEquals("Orthorhombic [111] vector incorrectly calculated", Math.sqrt(4+9+25), uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,1})), 1e-6);
		
		uc.updateCell(new Lattice.LatticeBuilder(8.28).setB(12.97).setC(7.15).setAl(91.05).setBe(116.26).setGa(90.15).build());
		assertEquals(15.21688, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,1})), 1e-5);
	}

}
