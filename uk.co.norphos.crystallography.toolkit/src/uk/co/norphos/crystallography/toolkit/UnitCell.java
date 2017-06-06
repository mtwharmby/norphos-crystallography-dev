package uk.co.norphos.crystallography.toolkit;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import uk.co.norphos.crystallography.api.IUnitCell;
import uk.co.norphos.crystallography.api.Lattice;

/**
 * FIXME
 * @author Michael Wharmby
 *
 */
public class UnitCell implements IUnitCell {
	
	private Lattice lattice;
	private IUnitCell reciprocal;
	private RealMatrix metricTensor, orthoMatrix, fracMatrix;
	private LUDecomposition metricTensorLUDecomp;

	/**
	 * FIXME
	 * @param realSpaceLattice
	 */
	public UnitCell(Lattice realSpaceLattice) {
		metricTensor = determineMetricTensor(realSpaceLattice);
		metricTensorLUDecomp = new LUDecomposition(metricTensor);
		lattice = CrystallographyFactory.createLattice(metricTensor);
		
		//Create the reciprocal space unit cell; the reciprocal of that is the present instance
		reciprocal = new UnitCell(metricTensorLUDecomp.getSolver().getInverse(), this);
		
		orthoMatrix = determineOrthogonalizationMatrix();
		LUDecomposition orthoMatLUDecomp = new LUDecomposition(orthoMatrix);
		fracMatrix = orthoMatLUDecomp.getSolver().getInverse();
	}
	
	/**
	 * FIXME
	 * @param metricTensor
	 */
	public UnitCell(RealMatrix metricTensor) {
		this(metricTensor, null);
	}
	
	private UnitCell(RealMatrix metricTensor, IUnitCell reciprocal) {
		this.metricTensor = metricTensor;
		metricTensorLUDecomp = new LUDecomposition(metricTensor);
		lattice = CrystallographyFactory.createLattice(metricTensor);
		
		if (reciprocal == null) {
			reciprocal = new UnitCell(metricTensorLUDecomp.getSolver().getInverse());
		} else {
			this.reciprocal = reciprocal;
		}
	}
	
	
	private RealMatrix determineMetricTensor(Lattice lattice) {
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
	
	private RealMatrix determineOrthogonalizationMatrix() {
		double[] elems = getConversionMatrixElements(lattice, getReciprocalLattice());
		return MatrixUtils.createRealMatrix(new double[][]{
			{elems[0], elems[1], elems[2]},
			{0       , elems[3], elems[4]},
			{0       , 0       , elems[5]}});
	}
	
	private double[] getConversionMatrixElements(Lattice one, Lattice two) {
		double[] elems = new double[6];
		elems[0] = one.getA();
		elems[1] = one.getB() * Math.cos(one.getGaR());
		elems[2] = one.getC() * Math.cos(one.getBeR());
		elems[3] = one.getB() * Math.sin(one.getGaR());
		elems[4] = -one.getC() * Math.sin(one.getBeR()) * Math.cos(two.getAlR());
		elems[5] = 1 / two.getC();
		return elems;
	}

	@Override
	public Lattice getLattice() {
		return lattice;
	}

	@Override
	public RealMatrix getMetricTensor() {
		return metricTensor;
	}

	@Override
	public IUnitCell getReciprocal() {
		return reciprocal;
	}
	
	@Override
	public RealMatrix getFractionalizationMatrix() {
		return fracMatrix;
	}

	@Override
	public RealMatrix getOrthogonalizationMatrix() {
		return orthoMatrix;
	}

	@Override
	public double calculateLength(Vector3D fracVec) {
		RealVector vector = new ArrayRealVector(fracVec.toArray());
		double product = vector.dotProduct(metricTensor.operate(vector));
		return Math.sqrt(product);
	}
	
	@Override
	public double calculateAngle(Vector3D fracVec1, Vector3D fracVec2) {
		RealVector vector1 = new ArrayRealVector(fracVec1.toArray());
		RealVector vector2 = new ArrayRealVector(fracVec2.toArray());
		
		double magVec1 = calculateLength(fracVec1);
		double magVec2 = calculateLength(fracVec2);
		
		return Math.acos(vector1.dotProduct(metricTensor.operate(vector2)) / (magVec1 * magVec2));
	}
	
	@Override
	public double calculateDihedralAngle(Vector3D site1, Vector3D site2, Vector3D site3, Vector3D site4) {
		Vector3D vector12 = site2.subtract(site1);
		Vector3D vector23 = site2.subtract(site3);
		Vector3D vector34 = site3.subtract(site4);
		
		Vector3D plane123 = latticeCrossProduct(vector12, vector23);
		Vector3D plane234 = latticeCrossProduct(vector34, vector23);
		return Math.acos(plane123.dotProduct(new Vector3D(metricTensor.operate(plane234.toArray()))) / (calculateLength(plane123) * calculateLength(plane234)));
		
		
	}
	
	public double latticeDotProduct(Vector3D vector1, Vector3D vector2) { //TODO add to API?
		Vector3D cartVec1 = orthogonalize(vector1);
		Vector3D cartVec2 = orthogonalize(vector2);
		return cartVec2.dotProduct(cartVec1);
	}
	
	public Vector3D latticeCrossProduct(Vector3D vector1, Vector3D vector2) { //TODO add to API?
		Vector3D cartVec1 = orthogonalize(vector1);
		Vector3D cartVec2 = orthogonalize(vector2);
		return fractionalize(cartVec1.crossProduct(cartVec2));
	}

	@Override
	public int compareTo(IUnitCell o) {
		// TODO Auto-generated method stub
		return 0;
	}
	
//	private Lattice determineReciprocalLattice() {
//		double rA = Math.sqrt(reciprocalMetricTensor.getEntry(0, 0));
//		double rB = Math.sqrt(reciprocalMetricTensor.getEntry(1, 1));
//		double rC = Math.sqrt(reciprocalMetricTensor.getEntry(2, 2));
//		double rAl = Math.toDegrees(Math.acos(reciprocalMetricTensor.getEntry(1, 2) / (rB * rC)));
//		double rBe = Math.toDegrees(Math.acos(reciprocalMetricTensor.getEntry(0, 2) / (rA * rC)));
//		double rGa = Math.toDegrees(Math.acos(reciprocalMetricTensor.getEntry(0, 1) / (rA * rB)));
//		
//		return new Lattice(rA, rB, rC, rAl, rBe, rGa, lattice.getPrincipleAxis());
//	}
//	

	
//	public Lattice getLattice() {
//		return lattice;
//	}
//	
//	public RealMatrix getMetricTensor() {
//		return metricTensor;
//	}
//	
//	public Lattice getReciprocalLattice() {
//		return reciprocalLattice;
//	}
//	
//	public RealMatrix getReciprocalMetricTensor() {
//		return reciprocalMetricTensor;
//	}
//	
//	public double getCellVolume() {
//		return volume;
//	}
//	
//	public double getReciprocalCellVolume() {
//		return reciprocalVolume;
//	}
//	
//	public double findVectorMagnitude(RealVector vector) {
//		return findVectorMagnitude(vector, metricTensor);
//	}
//	
//	public double findPlaneDSpacing(RealVector hklVector) {
//		return 1/findVectorMagnitude(hklVector, reciprocalMetricTensor);
//	}
//	
//	private double findVectorMagnitude(RealVector vector, RealMatrix tensor) {
//		double product = vector.dotProduct(tensor.operate(vector));
//		return Math.sqrt(product);
//	}

}
