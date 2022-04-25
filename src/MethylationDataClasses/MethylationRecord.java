package MethylationDataClasses;

import ngsep.genome.GenomicRegion;

public class MethylationRecord implements GenomicRegion {
	
	private String sequenceName;
	private int first;
	private int methylatedBaseCalls;
	private int unmethylatedBaseCalls;
	private int total;
	private double methCallPercentage;
	
	public MethylationRecord(String sequenceName, int first, int methylatedBaseCalls, int total) {
		this.sequenceName = sequenceName;
		this.first = first;
		this.methylatedBaseCalls = methylatedBaseCalls;
		this.total = total;
		this.unmethylatedBaseCalls = this.total - this.methylatedBaseCalls;
		this.methCallPercentage = (double) this.methylatedBaseCalls/this.total;
	}
	public boolean isCompatible(MethylationRecord methRecord) {
		return testCompatibility(this, methRecord);
	}
	@Override
	public String toString() {
		return "MethylationRecord [sequenceName=" + sequenceName + ", first=" + first + ", total=" + total + 
				 ", methylatedBaseCalls=" + methylatedBaseCalls + "]";
	}
	public boolean testCompatibility(MethylationRecord m1, MethylationRecord m2) {
		if(m1.equals(m2)) return true;
		if(m1.getFirst()!=m2.getFirst()) return false;
		if(m1.getSequenceName() != m2.getSequenceName()) return false;
		return true;
	}
	
	public String getSequenceName() {
		return sequenceName;
	}

	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}

	public int getFirst() {
		return first;
	}

	public void setFirst(int first) {
		this.first = first;
	}

	public int getMethylatedBaseCalls() {
		return methylatedBaseCalls;
	}

	public void setMethylatedBaseCalls(int methylatedBaseCalls) {
		this.methylatedBaseCalls = methylatedBaseCalls;
	}

	public int getUnmethylatedBaseCalls() {
		return unmethylatedBaseCalls;
	}

	public void setUnmethylatedBaseCalls(int unmethylatedBaseCalls) {
		this.unmethylatedBaseCalls = unmethylatedBaseCalls;
	}

	public int getTotal() {
		return total;
	}

	public void setTotal(int total) {
		this.total = total;
	}

	public double getMethCallPercentage() {
		return methCallPercentage;
	}

	public void setMethCallPercentage(double methCallPercentage) {
		this.methCallPercentage = methCallPercentage;
	}
	@Override
	public int getLast() {
		// TODO Auto-generated method stub
		return first + 1;
	}
	@Override
	public boolean isNegativeStrand() {
		// TODO Auto-generated method stub
		return false;
	}
	@Override
	public boolean isPositiveStrand() {
		// TODO Auto-generated method stub
		return false;
	}
	@Override
	public int length() {
		// TODO Auto-generated method stub
		return 1;
	}
}
