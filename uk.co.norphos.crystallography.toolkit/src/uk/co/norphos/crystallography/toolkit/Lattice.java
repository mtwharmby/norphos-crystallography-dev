package uk.co.norphos.crystallography.toolkit;

import java.util.ArrayList;
import java.util.List;

public class Lattice {

	private double a, b, c, al, be, ga;
	private double alR, beR, gaR; //Angles in radians for convenience
	private PrincipleAxis pAxis;

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
		this.pAxis = pAxis;
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
		return pAxis;
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

		public Lattice build() {
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
					b = c = a;
					if (angles[0] != 90) {
						//Rhombohedral
					} else {
						//Cubic
					}
				} else if ((a == b && a != c) || (b == null && c != null)) {
					//Tetragonal
					b = a;
					pAxis = PrincipleAxis.C;
				} else {
					//Orthorhombic (a != b != c)
					if (b == null || c == null) throw new IllegalArgumentException("Orthorhombic requires all three lengths to be given");
				}
			} else if (angles[2] == 120) {
				//Hexagonal
				b = a;
			} else if (!anglesCompared.contains(true)) {
				//Triclinic
				if (b == null || c == null) throw new IllegalArgumentException("Triclinic requires all three lengths to be given");
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

			//Principle axis is not meaningful for reciprocal lattice
			if (reciprocalLattice) pAxis = null;

			return new Lattice(a, b, c, angles[0], angles[1], angles[2], pAxis);
		}
	}

}
