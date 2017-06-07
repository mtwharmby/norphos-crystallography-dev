package uk.co.norphos.crystallography.api;

import java.io.Serializable;

public class MillerIndex implements Serializable, Comparable<MillerIndex> {

	private static final long serialVersionUID = 1162662794043478530L;

	private int[] indices;
	private double dSpacing;
	private String label;
	private double structureFactor, intensity;
	
	public MillerIndex(int h, int k, int l) {
		indices = new int[3];
		indices[0] = h;
		indices[1] = k;
		indices[2] = l;
	}
	
	public int[] getIndices() {
		return indices;
	}
	
	public int getH() {
		return indices[0];
	}
	
	public int getK() {
		return indices[1];
	}
	
	public int getL() {
		return indices[2];
	}
	
	public double getdSpacing() {
		return dSpacing;
	}

	public void setdSpacing(double dSpacing) {
		this.dSpacing = dSpacing;
	}

	public double getStructureFactor() {
		return structureFactor;
	}

	public void setStructureFactor(double structureFactor) {
		this.structureFactor = structureFactor;
	}

	public double getIntensity() {
		return intensity;
	}

	public void setIntensity(double intensity) {
		this.intensity = intensity;
	}

	@Override
	public int compareTo(MillerIndex arg0) {
		// TODO Auto-generated method stub
		//compare on d-spacing
		return 0;
	}

}
