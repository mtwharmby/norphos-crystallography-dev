package uk.co.norphos.crystallography.toolkit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class Lattice {

	private final double a, b, c, al, be, ga;
	private final double alR, beR, gaR; //Angles in radians for convenience
	private final PrincipleAxis principleAxis;

	private Lattice(Double a, Double b, Double c, Double al, Double be, Double ga, PrincipleAxis pAxis) {
		this.a = a;
		this.b = b;
		this.c = c;
		this.al = al;
		this.alR = Math.toRadians(al);
		this.be = be;
		this.beR = Math.toRadians(be);
		this.ga = ga;
		this.gaR = Math.toRadians(ga);
		this.principleAxis = pAxis;
	}

	public double a() {
		return a;
	}

	public double b() {
		return b;
	}

	public double c() {
		return c;
	}

	public double al() {
		return al;
	}

	public double be() {
		return be;
	}

	public double ga() {
		return ga;
	}

	public double alR() {
		return alR;
	}

	public double beR() {
		return beR;
	}

	public double gaR() {
		return gaR;
	}

	public PrincipleAxis getPrincipleAxis() {
		return principleAxis;
	}
	
	@Override
	public String toString() {
		return "Lattice [a=" + a + ", b=" + b + ", c=" + c + ", al=" + al 
				+ ", be=" + be + ", ga=" + ga + ", pAxis=" + principleAxis + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(a);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(al);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(alR);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(b);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(be);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(beR);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(c);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(ga);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(gaR);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((principleAxis == null) ? 0 : principleAxis.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Lattice other = (Lattice) obj;
		if (Double.doubleToLongBits(a) != Double.doubleToLongBits(other.a))
			return false;
		if (Double.doubleToLongBits(al) != Double.doubleToLongBits(other.al))
			return false;
		if (Double.doubleToLongBits(alR) != Double.doubleToLongBits(other.alR))
			return false;
		if (Double.doubleToLongBits(b) != Double.doubleToLongBits(other.b))
			return false;
		if (Double.doubleToLongBits(be) != Double.doubleToLongBits(other.be))
			return false;
		if (Double.doubleToLongBits(beR) != Double.doubleToLongBits(other.beR))
			return false;
		if (Double.doubleToLongBits(c) != Double.doubleToLongBits(other.c))
			return false;
		if (Double.doubleToLongBits(ga) != Double.doubleToLongBits(other.ga))
			return false;
		if (Double.doubleToLongBits(gaR) != Double.doubleToLongBits(other.gaR))
			return false;
		if (principleAxis != other.principleAxis)
			return false;
		return true;
	}



	public static class LatticeBuilder {

		private Double a, b, c;
		private Double[] angles = new Double[3];
		private boolean reciprocalLattice = false;
		private PrincipleAxis pAxis = PrincipleAxis.NONE;

		public LatticeBuilder(double a) {
			this.a = a;
		}

		public LatticeBuilder setB(double b) {
			this.b = b;
			return this;
		}

		public LatticeBuilder setC(double c) {
			this.c = c;
			return this;
		}

		public LatticeBuilder setAl(double al) {
			angles[0] = al;
			return this;
		}

		public LatticeBuilder setBe(double be) {
			angles[1] = be;
			return this;
		}

		public LatticeBuilder setGa(double ga) {
			angles[2] = ga;
			return this;
		}

		public LatticeBuilder setReciprocalLattice(boolean reciprocalLattice) {
			this.reciprocalLattice = reciprocalLattice;
			return this;
		}

		@Override
		public String toString() {
			return "LatticeBuilder [a=" + a + ", b=" + b + ", c=" + c + ", al=" + angles[0] 
					+ ", be=" + angles[1] + ", ga=" + angles[2] + ", pAxis=" + pAxis + "]";
		}

		public Lattice build() {
			if ((Collections.frequency(Arrays.asList(angles), null) == 2 && (b == null && c == null)) || 
					(angles[0] != null && angles[0] == angles[1] && angles[0] == angles[2])){
				//Rhombohedral
				if (b == null) b = a;
				if (c == null) c = a;
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
					if ((a == b && a == c) || (b == null && c == null)) {
						//Cubic
						b = c = a;
					} else if ((a == b && a != c) || (b == null && c != null)) {
						//Tetragonal
						b = a;
						pAxis = PrincipleAxis.C;
					} else {
						//Orthorhombic (a != b != c)
						if (b == null || c == null) throw new IllegalArgumentException("Orthorhombic requires all three lengths to be given");
					}
				} else if (!anglesCompared.contains(true)) {
					//Triclinic
					if (b == null || c == null) throw new IllegalArgumentException("Triclinic requires all three lengths to be given");
				} else {
					if (Arrays.asList(angles).contains(120.0) && (Collections.frequency(anglesCompared, true) == 1)) {
						//Hexagonal
						if (b == null) {
							b = a;
						}
						pAxis = PrincipleAxis.C;
					} else {
						//Monoclinic
						if (b == null || c == null) throw new IllegalArgumentException("Monoclinic requires all three lengths to be given");
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

			//Principle axis is not meaningful for reciprocal lattice
			if (reciprocalLattice) pAxis = null;

			return new Lattice(a, b, c, angles[0], angles[1], angles[2], pAxis);
		}
	}

}
