package uk.co.norphos.crystallography.toolkit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import uk.co.norphos.crystallography.api.CrystalSystem;
import uk.co.norphos.crystallography.api.Lattice;
import uk.co.norphos.crystallography.api.PrincipleAxis;

public class CrystallographyFactory {
	
	public static Lattice createLattice(double a, Double b, Double c, Double alpha, Double beta, Double gamma) {
		//TODO
		return null;
	}
	
	public static Lattice createLattice(RealMatrix metricTensor) {
		Double[] lengths = new Double[3], angles = new Double[3];
		int i,j,k;
		
		//Recover lengths...
		for (i = 0; i < 3; i++) {
			lengths[i] = Math.sqrt(metricTensor.getEntry(i, i));
		}
		//... and angles
		for (i = 0; i <3; i++) {
			j = (i + 1) % 3;
			k = (j + 1) % 3;
			
			angles[i] = Math.toDegrees(Math.acos(metricTensor.getEntry(j, k) / (lengths[j] * lengths[k])));
		}
		
		LUDecomposition metricTensorLUDecomp = new LUDecomposition(metricTensor);
		double volume = Math.sqrt(metricTensorLUDecomp.getDeterminant());
		
		return createLattice(lengths, angles, volume);
	}
	
	public static Lattice createLattice(Double[] lengths, Double[] angles, Double volume) {
		PrincipleAxis pAxis = PrincipleAxis.NONE;
		CrystalSystem crystalSystem;
		
		if ((Collections.frequency(Arrays.asList(angles), null) == 2 && (lengths[1] == null && lengths[2] == null)) || 
				(angles[0] != null && angles[0] == angles[1] && angles[0] == angles[2])){
			//Rhombohedral
			crystalSystem = CrystalSystem.RHOMBOHEDRAL;
			if (lengths[1] == null) lengths[1] = lengths[0];
			if (lengths[2] == null) lengths[2] = lengths[0];
			for (Double val : angles) {
				if (val != null) {
					angles[0] = val;
					angles[1] = val;
					angles[2] = val;
					break;
				}
			}
		} else {
			//Change all null angle values to 90 & compare them to one-another
			for (int i = 0; i < 3; i++) {
				if (angles[i] == null) angles[i] = 90.0;
			}
			List<Boolean> anglesCompared = new ArrayList<>(3); //[al-be, al-ga, be-ga]
			anglesCompared.add(angles[0].equals(angles[1]));
			anglesCompared.add(angles[0].equals(angles[2]));
			anglesCompared.add(angles[1].equals(angles[2]));

			if (!anglesCompared.contains(false)) {
				if ((lengths[0].equals(lengths[1]) && lengths[0].equals(lengths[2])) || (lengths[1] == null && lengths[2] == null)) {
					//Cubic
					crystalSystem = CrystalSystem.CUBIC;
					lengths[1] = lengths[2] = lengths[0];
				} else if ((lengths[0].equals(lengths[1]) && !lengths[0].equals(lengths[2])) || (lengths[1] == null && lengths[2] != null)) {
					//Tetragonal
					crystalSystem = CrystalSystem.TETRAGONAL;
					lengths[1] = lengths[0];
					pAxis = PrincipleAxis.C;
				} else {
					//Orthorhombic (a != b != c)
					crystalSystem = CrystalSystem.ORTHORHOMBIC;
					if (lengths[1] == null || lengths[2] == null) throw new IllegalArgumentException("Orthorhombic requires all three lengths to be given");
				}
			} else if (!anglesCompared.contains(true)) {
				//Triclinic
				crystalSystem = CrystalSystem.TRICLINIC;
				if (lengths[1] == null || lengths[2] == null) throw new IllegalArgumentException("Triclinic requires all three lengths to be given");
			} else {
				if (Arrays.asList(angles).contains(120.0) && (Collections.frequency(anglesCompared, true) == 1)) {
					//Hexagonal
					crystalSystem = CrystalSystem.HEXAGONAL;
					if (lengths[1] == null) {
						lengths[1] = lengths[0];
					}
					pAxis = PrincipleAxis.C;
				} else {
					//Monoclinic
					crystalSystem = CrystalSystem.MONOCLINIC;
					if (lengths[1] == null || lengths[2] == null) throw new IllegalArgumentException("Monoclinic requires all three lengths to be given");
					if (anglesCompared.get(0)) {
						//ga different
						pAxis = PrincipleAxis.C;
					} else if (anglesCompared.get(1)) {
						//be different
						pAxis = PrincipleAxis.B;
					} else {
						//al different
						pAxis = PrincipleAxis.A;
					}
				}
			}
		}
		return new Lattice(lengths, angles, volume, crystalSystem, pAxis);
				
	}
	



//	public static class LatticeBuilder {
//
//		private Double a, b, c;
//		private Double[] angles = new Double[3];
//		private boolean reciprocalLattice = false;
//		private PrincipleAxis pAxis = PrincipleAxis.NONE;
//
//		public LatticeBuilder(double a) {
//			this.a = a;
//		}
//
//		public LatticeBuilder setB(double b) {
//			this.b = b;
//			return this;
//		}
//
//		public LatticeBuilder setC(double c) {
//			this.c = c;
//			return this;
//		}
//
//		public LatticeBuilder setAl(double al) {
//			angles[0] = al;
//			return this;
//		}
//
//		public LatticeBuilder setBe(double be) {
//			angles[1] = be;
//			return this;
//		}
//
//		public LatticeBuilder setGa(double ga) {
//			angles[2] = ga;
//			return this;
//		}
//
//		public LatticeBuilder setReciprocalLattice(boolean reciprocalLattice) {
//			this.reciprocalLattice = reciprocalLattice;
//			return this;
//		}
//
//		@Override
//		public String toString() {
//			return "LatticeBuilder [a=" + a + ", b=" + b + ", c=" + c + ", al=" + angles[0] 
//					+ ", be=" + angles[1] + ", ga=" + angles[2] + ", pAxis=" + pAxis + "]";
//		}
//
//		public Lattice build() {
//			if ((Collections.frequency(Arrays.asList(angles), null) == 2 && (b == null && c == null)) || 
//					(angles[0] != null && angles[0] == angles[1] && angles[0] == angles[2])){
//				//Rhombohedral
//				if (b == null) b = a;
//				if (c == null) c = a;
//				for (Double val : angles) {
//					if (val != null) {
//						angles[0] = val;
//						angles[1] = val;
//						angles[2] = val;
//						break;
//					}
//				}
//			} else {
//				//Change all null angle values to 90 & compare them to one-another
//				for (int i = 0; i < 3; i++) {
//					if (angles[i] == null) angles[i] = 90.0;
//				}
//				List<Boolean> anglesCompared = new ArrayList<>(3); //[al-be, al-ga, be-ga]
//				anglesCompared.add(angles[0].equals(angles[1]));
//				anglesCompared.add(angles[0].equals(angles[2]));
//				anglesCompared.add(angles[1].equals(angles[2]));
//
//				if (!anglesCompared.contains(false)) {
//					if ((a == b && a == c) || (b == null && c == null)) {
//						//Cubic
//						b = c = a;
//					} else if ((a == b && a != c) || (b == null && c != null)) {
//						//Tetragonal
//						b = a;
//						pAxis = PrincipleAxis.C;
//					} else {
//						//Orthorhombic (a != b != c)
//						if (b == null || c == null) throw new IllegalArgumentException("Orthorhombic requires all three lengths to be given");
//					}
//				} else if (!anglesCompared.contains(true)) {
//					//Triclinic
//					if (b == null || c == null) throw new IllegalArgumentException("Triclinic requires all three lengths to be given");
//				} else {
//					if (Arrays.asList(angles).contains(120.0) && (Collections.frequency(anglesCompared, true) == 1)) {
//						//Hexagonal
//						if (b == null) {
//							b = a;
//						}
//						pAxis = PrincipleAxis.C;
//					} else {
//						//Monoclinic
//						if (b == null || c == null) throw new IllegalArgumentException("Monoclinic requires all three lengths to be given");
//						if (anglesCompared.get(0)) {
//							//ga different
//							pAxis = PrincipleAxis.C;
//						} else if (anglesCompared.get(1)) {
//							//be different
//							pAxis = PrincipleAxis.B;
//						} else {
//							//al different
//							pAxis = PrincipleAxis.A;
//						}
//					}
//				}
//			}
//
//			//Principle axis is not meaningful for reciprocal lattice
//			if (reciprocalLattice) pAxis = null;
//
//			return new Lattice(a, b, c, angles[0], angles[1], angles[2], pAxis);
//		}
//	}

}
