package MethylationDataClasses;

import java.util.SortedMap;

import ngsep.genome.GenomicRegion;

public class MethylatedRegion implements GenomicRegion {
	public static final String KEY_SEPARATOR = ",";
	
	private String sampleId;
	private String sequenceName;
	private int first;
	private int last;
	private int sampleCalls;
	private double methylatedControlPercentage;
	private double methylatedTreatmentPercentage;
	private double methylationPercentage;
	private int methylationCount;
	private double pValue;
	private double correctedPvalue;
	private boolean isDmr;

	private SortedMap<Integer, MethylationRecord[]> methylatedBases;
	
	//constructor for dmr
	public MethylatedRegion(String sequenceName, int first, int last,
			double methylatedControlPercentage, double methylatedTreatmentPercentage, double pValue, boolean isDmr,
			SortedMap<Integer, MethylationRecord[]> methylatedBases) {
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.isDmr = isDmr;
		this.sampleCalls = 0;
		this.methylatedControlPercentage = methylatedControlPercentage;
		this.methylatedTreatmentPercentage = methylatedTreatmentPercentage;
		this.methylationPercentage = methylatedTreatmentPercentage - methylatedControlPercentage;
		this.pValue = pValue;
		this.methylatedBases = methylatedBases;
	}
	
	//constructor for mr
	public MethylatedRegion(String sequenceName, int first, int last,
			double methylatedPercentage, double pValue,
			SortedMap<Integer, MethylationRecord[]> methylatedBases) {
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.sampleCalls = 0;
		this.isDmr = false;
		this.methylationPercentage = methylatedPercentage;
		this.pValue = pValue;
		this.methylatedBases = methylatedBases;
	}
	
	public MethylatedRegion(String sequenceName, int first, int last,
			int methylatedBaseCount, double pValue) {
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.sampleCalls = 0;
		this.isDmr = false;
		this.methylationCount = methylatedBaseCount;
		this.pValue = pValue;
	}
	
	public String decodeMRToString() {
		return this.getSequenceName() + KEY_SEPARATOR + this.getFirst() + KEY_SEPARATOR + this.getLast();
	}
	
	public boolean isDifferentiallyMethylated(){
		return isDmr;
		
	}
	
	@Override
	public int getFirst() {
		// TODO Auto-generated method stub
		return first;
	}

	@Override
	public int getLast() {
		// TODO Auto-generated method stub
		return last;
	}

	@Override
	public String getSequenceName() {
		// TODO Auto-generated method stub
		return sequenceName;
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
		return last - first + 1;
	}

	public double getMethylatedControlPercentage() {
		return methylatedControlPercentage;
	}

	public double getMethylatedTreatmentPercentage() {
		return methylatedTreatmentPercentage;
	}

	public double getMethylationPercentage() {
		return methylationPercentage;
	}

	public double getpValue() {
		return pValue;
	}

	public SortedMap<Integer, MethylationRecord[]> getMethylatedBases() {
		return methylatedBases;
	}
	public double getCorrectedPvalue() {
		return correctedPvalue;
	}
	
	private int getCalls() {
		// TODO Auto-generated method stub
		return sampleCalls;
	}

	public void setCorrectedPvalue(double correctedPvalue) {
		this.correctedPvalue = correctedPvalue;
	}
	public String getSampleId() {
		return sampleId;
	}

	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}
	public String toString() {
		if(isDmr) return this.getSequenceName() + "\t" + this.getFirst() + "\t" + this.getLast() + "\t" + this.getCalls() 
		+ "\t" + this.getMethylatedControlPercentage() + "\t" + this.getMethylatedTreatmentPercentage() +
				"\t" + this.getMethylationPercentage() + "\t" + this.getpValue() + "\t" + this.getCorrectedPvalue();
		else return this.getSequenceName() + "\t" + this.getFirst() + "\t" + this.getLast() + "\t" + this.getCalls() 
		+ "\t" + this.getMethylationPercentage() + "\t" + this.getpValue() + "\t" + this.getCorrectedPvalue();
	}
	
	public String simpleMRToString() {
		return this.getSequenceName() + "\t" + this.getFirst() + "\t" + this.getLast() + "\t" + this.getCalls() 
		+ "\t" + this.methylationCount;
		}

	public void setCalls(int calls) {
		// TODO Auto-generated method stub
		this.sampleCalls = calls;
	}
}
