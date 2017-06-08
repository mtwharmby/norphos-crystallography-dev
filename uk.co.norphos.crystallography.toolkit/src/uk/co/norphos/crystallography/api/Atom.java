package uk.co.norphos.crystallography.api;

import java.io.Serializable;
import java.util.Arrays;

public class Atom implements Serializable {

	private static final long serialVersionUID = 3959654678206070372L;

	private double[] coords;
	private String name, type;
	private double occ; 
	private Double radius;
	private double[][] uijMatrix;
	private Integer charge, isotope, coordinationNumber;

	public Atom(String name, String type, double x, double y, double z) {
		this(name, type, new double[]{x,y,z}, 1, new double[3][3], null, null, null, null);
		//FIXME uijMatrix should be equivalent to beq = 1 as default
	}

	public Atom(String name, String type, double x, double y, double z, double occ, double[][] uijMatrix) {
		this(name, type, new double[]{x,y,z}, occ, uijMatrix, null, null, null, null);
	}

	public Atom(String name, String type, double x, double y, double z, double occ, double[][] uijMatrix, 
			Double radius, Integer coordinationNumber, Integer charge, Integer isotope) {
		this(name, type, new double[]{x,y,z}, occ, uijMatrix, radius, coordinationNumber, charge, isotope);
	}

	public Atom(String name, String type, double[] coords, double occ, double[][] uijMatrix, 
			Double radius, Integer coordinationNumber, Integer charge, Integer isotope) {
		this.name = name;
		this.type = type;
		this.coords = coords;
		this.occ = occ;
		this.uijMatrix = uijMatrix;
		this.radius = radius;
		this.charge = charge;
		this.isotope = isotope;
	}

	public double[] getCoords() {
		return coords;
	}

	public void setCoords(double[] coords) {
		this.coords = coords;
	}

	public double getX() {
		return coords[0];
	}

	public void setX(double x) {
		coords[0] = x;
	}

	public double getY() {
		return coords[1];
	}

	public void setY(double y) {
		coords[1] = y;
	}

	public double getZ() {
		return coords[2];
	}

	public void setZ(double z) {
		coords[2] = z;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public double getOcc() {
		return occ;
	}

	public void setOcc(double occ) {
		this.occ = occ;
	}

	public double[][] getUijMatrix() {
		return uijMatrix;
	}

	public void setUijMatrix(double[][] uijMatrix) {
		this.uijMatrix = uijMatrix;
	}

	public double getRadius() {
		return radius;
	}

	public void setRadius(double radius) {
		this.radius = radius;
	}

	public Integer getCoordinationNumber() {
		return coordinationNumber;
	}

	public void setCoordinationNumber(Integer coordinationNumber) {
		this.coordinationNumber = coordinationNumber;
	}

	public int getCharge() {
		return charge;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public int getIsotope() {
		return isotope;
	}

	public void setIsotope(int isotope) {
		this.isotope = isotope;
	}

	@Override
	public String toString() {
		return "Atom [coords=" + Arrays.toString(coords) + ", name=" + name + ", type=" + type + ", occ=" + occ
				+ ", uijMatrix=" + Arrays.toString(uijMatrix) + "]";
		//FIXME Add radius, CN, charge & isotope through StringBuffer + if statement
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((charge == null) ? 0 : charge.hashCode());
		result = prime * result + ((coordinationNumber == null) ? 0 : coordinationNumber.hashCode());
		result = prime * result + Arrays.hashCode(coords);
		result = prime * result + ((isotope == null) ? 0 : isotope.hashCode());
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		long temp;
		temp = Double.doubleToLongBits(occ);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((radius == null) ? 0 : radius.hashCode());
		result = prime * result + ((type == null) ? 0 : type.hashCode());
		result = prime * result + Arrays.deepHashCode(uijMatrix);
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
		Atom other = (Atom) obj;
		if (charge == null) {
			if (other.charge != null)
				return false;
		} else if (!charge.equals(other.charge))
			return false;
		if (coordinationNumber == null) {
			if (other.coordinationNumber != null)
				return false;
		} else if (!coordinationNumber.equals(other.coordinationNumber))
			return false;
		if (!Arrays.equals(coords, other.coords))
			return false;
		if (isotope == null) {
			if (other.isotope != null)
				return false;
		} else if (!isotope.equals(other.isotope))
			return false;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (Double.doubleToLongBits(occ) != Double.doubleToLongBits(other.occ))
			return false;
		if (radius == null) {
			if (other.radius != null)
				return false;
		} else if (!radius.equals(other.radius))
			return false;
		if (type == null) {
			if (other.type != null)
				return false;
		} else if (!type.equals(other.type))
			return false;
		if (!Arrays.deepEquals(uijMatrix, other.uijMatrix))
			return false;
		return true;
	}
}
