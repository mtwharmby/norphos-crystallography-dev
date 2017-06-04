package uk.co.norphos.crystallography.toolkit;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.TestUtils;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Before;
import org.junit.Test;

import uk.co.norphos.crystallography.api.Lattice;
import uk.co.norphos.crystallography.api.LatticeException;

public class UnitCellTest {
	
	private UnitCell uc;
	
	private Lattice cubicLatt, orthoLatt, tricLatt, cubicRecipLatt, orthoRecipLatt, tricRecipLatt;
	private double cubicVol, orthoVol, tricVol;
	private RealMatrix cubicGMat, orthoGMat, tricGMat, cubicRecipGMat, orthoRecipGMat, tricRecipGMat;
	
	private Vector3D tricCoord1, tricCoord2, tricCoord3, tricCoord4;
	private double tricDistance12, tricDistance23, tricDistance34;
	private double tricAngle123, tricAngle234, tricDH1234;
	
	@Before
	public void setUp() {
		/*
		 * Structures from which data taken are located in the test-data directory
		 * Volumes calculated using VESTA
		 */
		cubicLatt =  new Lattice(5.43018, 5.43018, 5.43018, 90, 90, 90);
		cubicVol = 160.118936;
		cubicGMat = calculateMetricTensor(cubicLatt);
		cubicRecipLatt = calculateReciprocalLattice(cubicLatt, cubicVol);
		cubicRecipGMat = calculateMetricTensor(cubicRecipLatt);
		
		orthoLatt = new Lattice(23.49290, 6.34350, 19.63820, 90, 90, 90);
		orthoVol = 2926.626460;
		orthoGMat = calculateMetricTensor(orthoLatt);
		orthoRecipLatt = calculateReciprocalLattice(orthoLatt, orthoVol);
		orthoRecipGMat = calculateMetricTensor(orthoRecipLatt);
		
		tricLatt = new Lattice(7.19196, 8.12720, 8.12771, 82.4809, 69.2610, 69.2584);
		tricVol = 415.482298;
		tricGMat = calculateMetricTensor(tricLatt);
		tricRecipLatt = calculateReciprocalLattice(tricLatt, tricVol);
		tricRecipGMat = calculateMetricTensor(tricRecipLatt);
		tricCoord1 = new Vector3D(new double[]{0.61530, 0.02520, 0.07450});
		tricCoord2 = new Vector3D(new double[]{0.00810, 0.17050, 0.17120});
		tricDistance12 = 3.904;
		tricCoord3 = new Vector3D(0.25710, 0.95880, 0.60530);
		tricDistance23 = 8.673;
		tricAngle123 = 72.12;
		tricCoord4 = new Vector3D(0.38470, 0.97480, 0.92550);
		tricDistance34 = 3.082;
		tricAngle234 = 134.89;
		tricDH1234 = Math.toRadians(72.99);
	}
	
	private RealMatrix calculateMetricTensor(Lattice lattice) {
		double aa, ab, ac, bb, bc, cc;
		aa = Math.pow(lattice.getA(), 2);
		bb = Math.pow(lattice.getB(), 2);
		cc = Math.pow(lattice.getC(), 2);
		ab = lattice.getA() * lattice.getB() * Math.cos(lattice.getGaR());
		ac = lattice.getA() * lattice.getC() * Math.cos(lattice.getBeR());
		bc = lattice.getB() * lattice.getC() * Math.cos(lattice.getAlR());
		return MatrixUtils.createRealMatrix(new double[][]{{aa, ab, ac}, {ab, bb, bc}, {ac, bc, cc}});
	}
	
	private Lattice calculateReciprocalLattice(Lattice realLattice, double volume) {
		double aStar, bStar, cStar, alStar, beStar, gaStar;

		//Reciprocal lattice lengths
		aStar = realLattice.getB() * realLattice.getC() * Math.sin(realLattice.getAlR()) / volume;
		bStar = realLattice.getA() * realLattice.getC() * Math.sin(realLattice.getBeR()) / volume;
		cStar = realLattice.getA() * realLattice.getB() * Math.sin(realLattice.getGaR()) / volume;
		//Reciprocal lattice angles
		alStar = Math.toDegrees(Math.acos(
				(Math.cos(realLattice.getBeR()) * Math.cos(realLattice.getGaR()) - Math.cos(realLattice.getAlR())) / Math.abs(Math.sin(realLattice.getBeR()) * Math.sin(realLattice.getGaR()))
				));
		beStar = Math.toDegrees(Math.acos(
				(Math.cos(realLattice.getAlR()) * Math.cos(realLattice.getGaR()) - Math.cos(realLattice.getBeR())) / Math.abs(Math.sin(realLattice.getAlR()) * Math.sin(realLattice.getGaR()))
				));
		gaStar = Math.toDegrees(Math.acos(
				(Math.cos(realLattice.getAlR()) * Math.cos(realLattice.getBeR()) - Math.cos(realLattice.getGaR())) / Math.abs(Math.sin(realLattice.getAlR()) * Math.sin(realLattice.getBeR()))
				));
		return new Lattice(aStar, bStar, cStar, alStar, beStar, gaStar);
	}
	
	@Test
	public void testLatticeSetting() throws LatticeException {
		uc = new UnitCell(cubicLatt);
		assertEquals("Lattice is not the expected cubic", cubicLatt, uc.getLattice());
		TestUtils.assertEquals("Incorrect G-matrix for cubic", cubicGMat, uc.getMetricTensor(), 1e-10);
		assertEquals("Volume incorrect for cubic", cubicVol, uc.getVolume(), 5e-4);
		
		uc = new UnitCell(orthoLatt);
		assertEquals("Lattice is not the expected orthorhombic", orthoLatt, uc.getLattice());
		TestUtils.assertEquals("Incorrect G-matrix for orthorhombic", orthoGMat, uc.getMetricTensor(), 1e-10);
		assertEquals("Volume incorrect for orthorhombic", orthoVol, uc.getVolume(), 5e-4);
		
		uc = new UnitCell(tricLatt);
		assertEquals("Lattice is not the expected triclinic", tricLatt, uc.getLattice());
		TestUtils.assertEquals("Incorrect G-matrix", tricGMat, uc.getMetricTensor(), 1e-10);
		assertEquals("Volume incorrect", tricVol, uc.getVolume(), 5e-4);
	}
	
	@Test
	public void testReciprocalLatticeCalcs() {
		uc = new UnitCell(cubicLatt);
		TestUtils.assertEquals("Incorrect reciprocal metric tensor for cubic", cubicRecipGMat, uc.getReciprocalMetricTensor(), 5e-9);
		
		uc = new UnitCell(orthoLatt);
		TestUtils.assertEquals("Incorrect reciprocal metric tensor for orthorhombic", orthoRecipGMat, uc.getReciprocalMetricTensor(), 5e-9);
		
		uc = new UnitCell(tricLatt);
		TestUtils.assertEquals("Incorrect reciprocal metric tensor for triclinic", tricRecipGMat, uc.getReciprocalMetricTensor(), 5e-9);
	}
	
	@Test
	public void testDistanceCalculation() {
		uc = new UnitCell(tricLatt);
				
		Vector3D fracVec1 = tricCoord2.subtract(tricCoord1);
		assertEquals("Wrong distance between sites 1 & 2 by vector", tricDistance12, uc.calculateLength(fracVec1), 5e-4);
		assertEquals("Wrong distance between sites 2 & 3 by sites", tricDistance23, uc.calculateDistance(tricCoord2, tricCoord3), 5e-4);
	}
	
	@Test
	public void testAngleCalculation() {
		uc = new UnitCell(tricLatt);
		
		Vector3D fracVec1, fracVec2;
		fracVec1 = new Vector3D(new double[]{1,0,0});
		fracVec2 = new Vector3D(new double[]{0,0,1});
		assertEquals("Expecting beta angle", tricLatt.getBeR(), uc.calculateAngle(fracVec1, fracVec2), 1e-4);
		
		fracVec1 = tricCoord2.subtract(tricCoord1);
		fracVec2 = tricCoord2.subtract(tricCoord3);
		assertEquals("Wrong angle between sites 1/2/3 by vector", tricAngle123, Math.toDegrees(uc.calculateAngle(fracVec1, fracVec2)), 5e-3);
		assertEquals("Wrong angle between sites 2/3/4 by sites", tricAngle234, Math.toDegrees(uc.calculateAngle(tricCoord2, tricCoord3, tricCoord4)), 5e-3);
	}
	
	
	
	
	
	
//	@Test
//	public void testVectorMagnitude() {
//		uc = new UnitCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).build());
//		assertEquals("Orthorhombic [100] vector incorrectly calculated", 2.0, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,0,0})), 1e-6);
//		assertEquals("Orthorhombic [010] vector incorrectly calculated", 3.0, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{0,1,0})), 1e-6);
//		assertEquals("Orthorhombic [110] vector incorrectly calculated", Math.sqrt(4+9), uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,0})), 1e-6);
//		assertEquals("Orthorhombic [111] vector incorrectly calculated", Math.sqrt(4+9+25), uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,1})), 1e-6);
//		
//		uc.updateCell(new Lattice.LatticeBuilder(8.28).setB(12.97).setC(7.15).setAl(91.05).setBe(116.26).setGa(90.15).build());
//		assertEquals(15.21688, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,1})), 1e-5);
//	}
//	
//	@Test
//	public void testLatticePlaneDSpacing() {
//		//Try more complicated triclinic case: This is Anorthoclase
//		Lattice anorthoclase = new Lattice.LatticeBuilder(8.28).setB(12.97).setC(7.15).setAl(91.05).setBe(116.26).setGa(90.15).build();
//		uc = new UnitCell(anorthoclase);
//		
//		//d-spacing of planes (calculated with PowderCell)
//		//(100) 7.42494; (010) 12.9669; (001) 6.41057; (110) 6.41039; (111) 3.84084
//		//Could do this on the fly too...
//		assertEquals("(100) spacing incorrect", 7.42494, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{1,0,0})), 1e-5);
//		assertEquals("(010) spacing incorrect", 12.96689, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{0,1,0})), 1e-5);
//		assertEquals("(001) spacing incorrect", 6.41057, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{0,0,1})), 1e-5);
//		assertEquals("(110) spacing incorrect", 6.41039, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{1,1,0})), 1e-5);
//		assertEquals("(111) spacing incorrect", 3.84084, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{1,1,1})), 1e-5);
//	}
	
	
	
	
	
	
	
	
	
//	@Test
//	public void testMetricTensor() {
//		//Cubic a = 3
//		uc = new UnitCell(new Lattice.LatticeBuilder(3).build());
//		RealMatrix fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{9,0,0},{0,9,0},{0,0,9}});
//		TestUtils.assertEquals("Cubic tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
//		
//		//Rhombohedral a = 3; al = 60
//		uc.updateCell(new Lattice.LatticeBuilder(3).setAl(60).build());
//		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{9,4.5,4.5},{4.5,9,4.5},{4.5,4.5,9}});
//		TestUtils.assertEquals("Rhombohedral tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
//
//		//Hexagonal a = 5; c = 2; ga = 120
//		uc.updateCell(new Lattice.LatticeBuilder(5).setC(2).setGa(120).build());
//		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{25,-12.5,0},{-12.5,25,0},{0,0,4}});
//		TestUtils.assertEquals("Hexagonal tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
//		
//		//Tetragonal a = b = 2; c = 5
//		uc.updateCell(new Lattice.LatticeBuilder(2).setC(5).build());
//		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{4,0,0},{0,4,0},{0,0,25}});
//		TestUtils.assertEquals("Tetragonal tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
//		
//		//Orthorhombic a = 5; b = 3; c = 2
//		uc.updateCell(new Lattice.LatticeBuilder(5).setB(3).setC(2).build());
//		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{25,0,0},{0,9,0},{0,0,4}});
//		TestUtils.assertEquals("Orthorhombic  tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
//		
//		//Monoclinic a = 2; b = 3; c = 5; be = 30deg
//		uc.updateCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).setBe(30).build());
//		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{4,0,5*Math.sqrt(3)},{0,9,0},{5*Math.sqrt(3),0,25}});
//		TestUtils.assertEquals("Monoclinic tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
//		
//		//Triclinic #a = 3; b = 5; c = 2; al = 30deg; be = 45 deg; ga = 60 deg
//		uc.updateCell(new Lattice.LatticeBuilder(3).setB(5).setC(2).setAl(30).setBe(45).setGa(60).build());
//		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{9,7.5,6/Math.sqrt(2)},{7.5, 25, 5*Math.sqrt(3)},{6/Math.sqrt(2), 5*Math.sqrt(3), 4}});
//		TestUtils.assertEquals("Triclinic tensor incorrectly calculated", fakeTensor, uc.getMetricTensor(), 1e-10);
//	}
//	
//	@Test
//	public void testCellVolume() {
//		uc = new UnitCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).build());
//		assertEquals(30.0, uc.getCellVolume(), 0.05);
//	}
//	
//	@Test
//	public void testReciprocalMetricTensor() {
//		uc = new UnitCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).build());
//		RealMatrix fakeTensor = MatrixUtils.createRealMatrix(new double[][]{{4,0,0}, {0,9,0}, {0,0,25}});
//		TestUtils.assertEquals("Orthorhombic tensor incorrectly calculated",  fakeTensor, uc.getMetricTensor(), 1e-10);
//		
//		RealMatrix fakeReciprocalTensor = MatrixUtils.createRealMatrix(new double[][]{
//			{9.*25./Math.pow(uc.getCellVolume(),2), 0, 0},
//			{0, 4.*25./Math.pow(uc.getCellVolume(),2), 0},
//			{0, 0, 4.*9./Math.pow(uc.getCellVolume(), 2)}});
//		TestUtils.assertEquals("Orthorhombic reciprocal tensor incorrectly calculated", fakeReciprocalTensor, uc.getReciprocalMetricTensor(), 1e-10);
//		
//		//Try more complicated triclinic case: This is Anorthoclase
//		Lattice anorthoclase = new Lattice.LatticeBuilder(8.28).setB(12.97).setC(7.15).setAl(91.05).setBe(116.26).setGa(90.15).build();
//		uc.updateCell(anorthoclase);
//		//Real space metric tensor
//		double g00 = Math.pow(anorthoclase.a(),2);
//		double g11 = Math.pow(anorthoclase.b(),2);
//		double g22 = Math.pow(anorthoclase.c(),2);
//		double g01 = anorthoclase.a()*anorthoclase.b()*Math.cos(anorthoclase.gaR());
//		double g02 = anorthoclase.a()*anorthoclase.c()*Math.cos(anorthoclase.beR());
//		double g12 = anorthoclase.b()*anorthoclase.c()*Math.cos(anorthoclase.alR());
//		fakeTensor = MatrixUtils.createRealMatrix(new double[][]{
//			{g00, g01, g02},
//			{g01, g11, g12},
//			{g02, g12, g22}});
//		TestUtils.assertEquals("Anorthoclase tensor incorrectly calculated",  fakeTensor, uc.getMetricTensor(), 1e-10);
//		
//		//Reciprocal space metric tensor
//		g00 = Math.pow(anorthoclase.b(),2) * Math.pow(anorthoclase.c(),2) * Math.pow(Math.sin(anorthoclase.alR()),2) / Math.pow(uc.getCellVolume(), 2);
//		g11 = Math.pow(anorthoclase.a(),2) * Math.pow(anorthoclase.c(),2) * Math.pow(Math.sin(anorthoclase.beR()),2) / Math.pow(uc.getCellVolume(), 2);
//		g22 = Math.pow(anorthoclase.a(),2) * Math.pow(anorthoclase.b(),2) * Math.pow(Math.sin(anorthoclase.gaR()),2) / Math.pow(uc.getCellVolume(), 2);
//		g01 = anorthoclase.a() * anorthoclase.b() * Math.pow(anorthoclase.c(),2) * (Math.cos(anorthoclase.alR()) * Math.cos(anorthoclase.beR()) - Math.cos(anorthoclase.gaR())) / Math.pow(uc.getCellVolume(),2);
//		g02 = anorthoclase.a() * Math.pow(anorthoclase.b(),2) * anorthoclase.c() * (Math.cos(anorthoclase.alR()) * Math.cos(anorthoclase.gaR()) - Math.cos(anorthoclase.beR())) / Math.pow(uc.getCellVolume(),2);
//		g12 = Math.pow(anorthoclase.a(),2) * anorthoclase.b() * anorthoclase.c() * (Math.cos(anorthoclase.beR()) * Math.cos(anorthoclase.gaR()) - Math.cos(anorthoclase.alR())) / Math.pow(uc.getCellVolume(),2);
//		fakeReciprocalTensor =  MatrixUtils.createRealMatrix(new double[][]{
//			{g00, g01, g02},
//			{g01, g11, g12},
//			{g02, g12, g22}});
//		TestUtils.assertEquals("Anorthoclase reciprocal tensor incorrectly calculated",  fakeReciprocalTensor, uc.getReciprocalMetricTensor(), 1e-10);
//	}
//	
//	@Test
//	public void testVectorMagnitude() {
//		uc = new UnitCell(new Lattice.LatticeBuilder(2).setB(3).setC(5).build());
//		assertEquals("Orthorhombic [100] vector incorrectly calculated", 2.0, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,0,0})), 1e-6);
//		assertEquals("Orthorhombic [010] vector incorrectly calculated", 3.0, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{0,1,0})), 1e-6);
//		assertEquals("Orthorhombic [110] vector incorrectly calculated", Math.sqrt(4+9), uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,0})), 1e-6);
//		assertEquals("Orthorhombic [111] vector incorrectly calculated", Math.sqrt(4+9+25), uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,1})), 1e-6);
//		
//		uc.updateCell(new Lattice.LatticeBuilder(8.28).setB(12.97).setC(7.15).setAl(91.05).setBe(116.26).setGa(90.15).build());
//		assertEquals(15.21688, uc.findVectorMagnitude(MatrixUtils.createRealVector(new double[]{1,1,1})), 1e-5);
//	}
//	
//	@Test
//	public void testLatticePlaneDSpacing() {
//		//Try more complicated triclinic case: This is Anorthoclase
//		Lattice anorthoclase = new Lattice.LatticeBuilder(8.28).setB(12.97).setC(7.15).setAl(91.05).setBe(116.26).setGa(90.15).build();
//		uc = new UnitCell(anorthoclase);
//		
//		//d-spacing of planes (calculated with PowderCell)
//		//(100) 7.42494; (010) 12.9669; (001) 6.41057; (110) 6.41039; (111) 3.84084
//		//Could do this on the fly too...
//		assertEquals("(100) spacing incorrect", 7.42494, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{1,0,0})), 1e-5);
//		assertEquals("(010) spacing incorrect", 12.96689, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{0,1,0})), 1e-5);
//		assertEquals("(001) spacing incorrect", 6.41057, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{0,0,1})), 1e-5);
//		assertEquals("(110) spacing incorrect", 6.41039, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{1,1,0})), 1e-5);
//		assertEquals("(111) spacing incorrect", 3.84084, uc.findPlaneDSpacing(MatrixUtils.createRealVector(new double[]{1,1,1})), 1e-5);
//	}
//		
//	@Test
//	public void testGetReciprocalCell() {
//		Lattice orthoLat = new Lattice.LatticeBuilder(2).setB(3).setC(5).build();
//		uc = new UnitCell(orthoLat);
//		Lattice fakeReciprocalLat = calculateReciprocalLattice(orthoLat, uc.getCellVolume());
//		//Compare lattices
//		assertEquals("a* incorrectly calculated", fakeReciprocalLat.a(), uc.getReciprocalLattice().a(), 1e-10);
//		assertEquals("b* incorrectly calculated", fakeReciprocalLat.b(), uc.getReciprocalLattice().b(), 1e-10);
//		assertEquals("c* incorrectly calculated", fakeReciprocalLat.c(), uc.getReciprocalLattice().c(), 1e-10);
//		assertEquals("al* incorrectly calculated", fakeReciprocalLat.al(), uc.getReciprocalLattice().al(), 1e-10);
//		assertEquals("be* incorrectly calculated", fakeReciprocalLat.be(), uc.getReciprocalLattice().be(), 1e-10);
//		assertEquals("ga* incorrectly calculated", fakeReciprocalLat.ga(), uc.getReciprocalLattice().ga(), 1e-10);
//		
//		Lattice anorthoclase = new Lattice.LatticeBuilder(8.28).setB(12.97).setC(7.15).setAl(91.05).setBe(116.26).setGa(90.15).build();
//		uc.updateCell(anorthoclase);
//		fakeReciprocalLat = calculateReciprocalLattice(anorthoclase, uc.getCellVolume());
//		assertEquals("a* incorrectly calculated", fakeReciprocalLat.a(), uc.getReciprocalLattice().a(), 1e-10);
//		assertEquals("b* incorrectly calculated", fakeReciprocalLat.b(), uc.getReciprocalLattice().b(), 1e-10);
//		assertEquals("c* incorrectly calculated", fakeReciprocalLat.c(), uc.getReciprocalLattice().c(), 1e-10);
//		assertEquals("al* incorrectly calculated", fakeReciprocalLat.al(), uc.getReciprocalLattice().al(), 1e-10);
//		assertEquals("be* incorrectly calculated", fakeReciprocalLat.be(), uc.getReciprocalLattice().be(), 1e-10);
//		assertEquals("ga* incorrectly calculated", fakeReciprocalLat.ga(), uc.getReciprocalLattice().ga(), 1e-10);
//	}
//	
//	private Lattice calculateReciprocalLattice(Lattice realLattice, double volume) {
//		double aStar, bStar, cStar, alStar, beStar, gaStar;
//		
//		//Reciprocal lattice lengths
//		aStar = realLattice.b() * realLattice.c() * Math.sin(realLattice.alR()) / volume;
//		bStar = realLattice.a() * realLattice.c() * Math.sin(realLattice.beR()) / volume;
//		cStar = realLattice.a() * realLattice.b() * Math.sin(realLattice.gaR()) / volume;
//		//Reciprocal lattice angles
//		alStar = Math.toDegrees(Math.acos(
//				(Math.cos(realLattice.beR()) * Math.cos(realLattice.gaR()) - Math.cos(realLattice.alR())) / Math.abs(Math.sin(realLattice.beR()) * Math.sin(realLattice.gaR()))
//				));
//		beStar = Math.toDegrees(Math.acos(
//				(Math.cos(realLattice.alR()) * Math.cos(realLattice.gaR()) - Math.cos(realLattice.beR())) / Math.abs(Math.sin(realLattice.alR()) * Math.sin(realLattice.gaR()))
//				));
//		gaStar = Math.toDegrees(Math.acos(
//				(Math.cos(realLattice.alR()) * Math.cos(realLattice.beR()) - Math.cos(realLattice.gaR())) / Math.abs(Math.sin(realLattice.alR()) * Math.sin(realLattice.beR()))
//						));
//		return new Lattice.LatticeBuilder(aStar).setB(bStar).setC(cStar).setAl(alStar).setBe(beStar).setGa(gaStar).build();
//	}
}