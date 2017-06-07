package uk.co.norphos.crystallography.api;

import java.io.Serializable;

public class Atom implements Serializable {

	private static final long serialVersionUID = 3959654678206070372L;
	
	private double[] coords;
	private String name, type;
	private double occ, radius;
	private double[][] uijMatrix;
	private int charge, isotope;

	public Atom(double x, double y, double z) {
		coords = new double[3];
		coords[0] = x;
		coords[1] = y;
		coords[2] = z;
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

	public double getRadius() {
		return radius;
	}

	public void setRadius(double radius) {
		this.radius = radius;
	}

	public double[][] getUijMatrix() {
		return uijMatrix;
	}

	public void setUijMatrix(double[][] uijMatrix) {
		this.uijMatrix = uijMatrix;
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
}
