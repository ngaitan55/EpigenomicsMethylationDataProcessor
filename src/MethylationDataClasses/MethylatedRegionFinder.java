package MethylationDataClasses;

import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.BetaDistribution;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.math.Distribution;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class MethylatedRegionFinder {
	
	public static final double METHYLATION_POS_TRESHOLD = 0.3;
	public static final double BETA_DISTRIBUTION_ALPHA_PARAMETER = 0.5;
	public static final double BETA_DISTRIBUTION_BETA_PARAMETER = 30;
	
	private ReferenceGenome referenceGenome;
	private double significanceLevel;
	private BetaDistribution betaDistribution;
	private int slidingWindowSize;
	private String [] customSampleIds;
	private PrintWriter writer = null;
	private QualifiedSequenceList sequences;
	private List<Distribution> methylationPercentageDistributions;
	private GenomicRegionSortedCollection<MethylatedRegion> methylatedRegions;
	private Map<String, Integer> mrSampleCounts; 
	
	public MethylatedRegionFinder(String genomeFile, int swLength, double initialAlpha) throws IOException {
		mrSampleCounts = new HashMap<>();
		significanceLevel = initialAlpha;
		slidingWindowSize = swLength;
		referenceGenome = new ReferenceGenome(genomeFile);
		sequences = referenceGenome.getSequencesList();
		betaDistribution = new BetaDistribution(BETA_DISTRIBUTION_ALPHA_PARAMETER, BETA_DISTRIBUTION_BETA_PARAMETER);
		methylationPercentageDistributions = new ArrayList<>();
		methylatedRegions = new GenomicRegionSortedCollection<>(sequences);
	}
	
	public void run(List<String> files, String outPrefix) throws IOException {
		writer = new PrintWriter(outPrefix + ".MR.tsv");
		printHeaderFields();
		//readAndProcessDistributions(files);
		readAndProcessSamples(files);
		//System.out.println("Applying correction");
		//System.out.println("Printing output file");
		printMrs();
	}
	
	public void readAndProcessDistributions(List<String> files) throws IOException{
		List<MethylationSampleFileReader> readers = new ArrayList<>();
		List<Iterator<MethylationRecord>> iterators = new ArrayList<>();
		QualifiedSequenceList sequences = referenceGenome.getSequencesList();
		for (String file:files) {
			MethylationSampleFileReader reader = new MethylationSampleFileReader(file);
			readers.add(reader);
			iterators.add(reader.iterator());
		}
		int nSamples = files.size();
		int currentSWBegin = 0;
		int currentSWLimit = slidingWindowSize;
		for(int r = 0; r < nSamples; r++){
			//MethylationRecord[] currentNext = new MethylationRecord[nSamples];
			Iterator<MethylationRecord> it = iterators.get(r);
			Distribution sampleDistribution = new Distribution(0, 1, 0.01);
			for(QualifiedSequence seq : sequences){
				int seqSize = seq.getLength();
				int nWindows = seqSize/slidingWindowSize;
				int lastWindowSize = seqSize%slidingWindowSize;
				String seqName = seq.getName();
				boolean iteratorOnHold = false;
				MethylationRecord nextRecord = null;
				for(int i = 0; i < nWindows; i++) {
					currentSWLimit = i + 1 == nWindows ? lastWindowSize : currentSWLimit;
					GenomicRegion slidingWindow = new GenomicRegionImpl(seqName, currentSWBegin, currentSWLimit);
					int windowMethylationCount = 0;
					while(it.hasNext()) {
						if(!iteratorOnHold) nextRecord = it.next();
						int cmp = compare(slidingWindow, nextRecord);
						if(cmp == -1) {
							iteratorOnHold = true;
							break;
						}
						else if(cmp == 1) {
							iteratorOnHold = false;
							if(!it.hasNext()) 
								break;
						}
						else {
							int methylated = nextRecord.getMethylatedBaseCalls();
							int total = nextRecord.getTotal();
							if(baseIsMethylated(methylated, total)) {
								windowMethylationCount++;
							}
							iteratorOnHold = false;
						}
					}
					double windowMethylationPercentage = (double) windowMethylationCount/slidingWindowSize;
					currentSWBegin += slidingWindowSize;
					currentSWLimit += slidingWindowSize;
					sampleDistribution.processDatapoint( (double) windowMethylationPercentage);
				}
			}
			currentSWBegin = 0;
			currentSWLimit = slidingWindowSize;
			methylationPercentageDistributions.add(sampleDistribution);
		}
	}
	
	public void readAndProcessSamples(List<String> files) throws IOException {
		List<MethylationSampleFileReader> readers = new ArrayList<>();
		List<Iterator<MethylationRecord>> iterators = new ArrayList<>();
		QualifiedSequenceList sequences = referenceGenome.getSequencesList();
		for (String file:files) {
			MethylationSampleFileReader reader = new MethylationSampleFileReader(file);
			readers.add(reader);
			iterators.add(reader.iterator());
		}
		int nSamples = files.size();
		int currentSWBegin = 0;
		int currentSWLimit = slidingWindowSize;
		for(int r = 0; r < nSamples; r++){
			//MethylationRecord[] currentNext = new MethylationRecord[nSamples];
			int nHypotheses = 0;
			int nWindows = 0;
			List<MethylatedRegion> sampleMrs = new ArrayList<>();
			Iterator<MethylationRecord> it = iterators.get(r);
			//Distribution sampleDistribution = methylationPercentageDistributions.get(r);
			for(QualifiedSequence seq : sequences){
				int seqSize = seq.getLength();
				nWindows = seqSize/slidingWindowSize;
				nHypotheses += nWindows;
				int lastWindowSize = seqSize%slidingWindowSize;
				String seqName = seq.getName();
				//System.out.println("Processing sequence: " + seqName + " nW=" + nWindows);
				boolean iteratorOnHold = false;
				MethylationRecord nextRecord = null;
				for(int i = 0; i < nWindows; i++) {
					currentSWLimit = i + 1 == nWindows ? lastWindowSize : currentSWLimit;
					GenomicRegion slidingWindow = new GenomicRegionImpl(seqName, currentSWBegin, currentSWLimit);
					int windowMethylationCount = 0;
					//SortedMap<Integer, MethylationRecord[]> dmrRecords = new TreeMap<>();
					while(it.hasNext()) {
						if(!iteratorOnHold) nextRecord = it.next();
						int cmp = compare(slidingWindow, nextRecord);
						//System.out.println("cmp=" + cmp);
						if(cmp == -1) {
							iteratorOnHold = true;
							break;
						}
						else if(cmp == 1) {
							iteratorOnHold = false;
							//if(!it.hasNext()) 
							break;
						}
						else {
							int methylated = nextRecord.getMethylatedBaseCalls();
							int total = nextRecord.getTotal();
							if(baseIsMethylated(methylated, total)) {
								windowMethylationCount++;
							}
							iteratorOnHold = false;
						}
					}
					//while(seqName == nextRecord.getSequenceName()) {
						//if(it.hasNext()) nextRecord = it.next();
					//}
					double windowMethylationPercentage = (double) windowMethylationCount/slidingWindowSize;
					//double pValue = sampleDistribution.getEmpiricalPvalue( (double) windowMethylationPercentage);
					double pValue = 1 - betaDistribution.cumulativeProbability(windowMethylationPercentage);
					//if(r == 0 || r == 6 || r == 10 || r == 15) System.out.println(pValue);
					if(testWindowMethylation(pValue)) {
						MethylatedRegion mr = new MethylatedRegion(seqName, slidingWindow.getFirst(), 
								slidingWindow.getLast(), windowMethylationPercentage*100, pValue, null);
						sampleMrs.add(mr);
						//System.out.println("mrFirst=" + mr.getFirst() + " CHR=" + mr.getSequenceName());
						String mrKey = mr.decodeMRToString();
						mrSampleCounts.compute(mrKey, (k,v) -> (v==null) ? 1:v+1);
					}
					currentSWBegin += slidingWindowSize;
					currentSWLimit += slidingWindowSize;
				}
				//System.out.println("nPerChr=" + sampleMrs.size() + " CHR=" + seqName);
				currentSWBegin = 0;
				currentSWLimit = slidingWindowSize;
			}
			//System.out.println("nh="+nHypotheses[r]++);
			applyBenjaminiHochbergCorrection(sampleMrs, nHypotheses);
			//applyBonferroniCorrection(sampleMrs, nHypotheses);
			methylatedRegions.addAll(sampleMrs);
			//Print the methylation bases distributions
			//String currentDistFile = "distribution" + r + ".txt";
			//PrintStream distPrinter = new PrintStream(currentDistFile);
			//sampleDistribution.printDistribution(new PrintStream(distPrinter));
			//distPrinter.close();
			if(r == 0 || r == 6 || r == 10 || r == 15) System.out.println("#");
			//System.out.println("Sample=" + r);
		}
		//List<String> chrs =  methylatedRegions.getSequenceNames().getNamesStringList();
		//for(String chr : chrs) {
			//System.out.println("CHR=" + chr);
		//}
	}
	
	public void applyBonferroniCorrection (List<MethylatedRegion> sampleMrs, int nHypotheses) {
		//double m = (double) dmrs.size();
		double m = (double) nHypotheses;
		double alpha =  significanceLevel / m;
		//System.out.println("m=" + m + " alpha_corrected=" + alpha);
		List<MethylatedRegion> toRemove = new ArrayList<>();
		for (MethylatedRegion mr : sampleMrs) {
			if(mr.getpValue() > alpha) toRemove.add(mr);
			else mr.setCorrectedPvalue(mr.getpValue() * m);
		}
		sampleMrs.removeAll(toRemove);
	}
	
	private void applyBenjaminiHochbergCorrection(List<MethylatedRegion> sampleMrs, int nHypotheses) {
		// TODO Auto-generated method stub
		Collections.sort(sampleMrs, Comparator.comparingDouble(mr -> mr.getpValue()));
		int m = nHypotheses;
		double fdr = significanceLevel;
		List<MethylatedRegion> toKeep = new ArrayList<>();
		int k = 0;
		while(k < sampleMrs.size()) {
			MethylatedRegion mr = sampleMrs.get(k);
			double pValue = mr.getpValue();
			double currentCriticalValue = calculateCorrectedBHCriticalValue(k + 1, m, fdr);
			//System.out.println("$pvalue=" + pValue + " alpha=" + currentCriticalValue);
			if(pValue > currentCriticalValue) {
				//System.out.println("bool " + (pValue > currentCriticalValue));
				break;
			}
			else {
				mr.setCorrectedPvalue(currentCriticalValue);
				toKeep.add(mr);
			}
			k++;
		}
		sampleMrs.retainAll(toKeep);
	}
	
	private static double calculateCorrectedBHCriticalValue(int i, int n, double fdr) {
		double alpha = (double) i/n;
		//System.out.println("i=" + i + " n=" + n + " alpha=" + alpha);
		//System.out.println("bjalpha=" + alpha*fdr);
		return alpha*fdr;
	}
	
	private boolean testWindowMethylation(double pValue) {
		// TODO Auto-generated method stub
		return pValue < significanceLevel;
	}

	public int compare(GenomicRegion sw, MethylationRecord m) {
		int p1 = sequences.indexOf(sw.getSequenceName());
		int p2 = sequences.indexOf(m.getSequenceName());
		//System.out.println("sw=" + sw.getSequenceName() + " " + sw.getFirst() + " m=" + m.getSequenceName() + " " + m.getFirst());
		if(p1 < p2) return -1;
		//This case should not happen
		if (p1>p2) return 1;
		if(m.getFirst() > sw.getLast()) return -1;
		if(m.getFirst() >=  sw.getFirst() && m.getFirst() <= sw.getLast()) return 0;
		//This case should not happen
		if(sw.getFirst() > m.getFirst()) return 1;
		return 0;
	}
	
	public boolean baseIsMethylated(int methylated, int total) {
		return (double) methylated / total >= METHYLATION_POS_TRESHOLD; 
	}
	
	public void printHeaderFields() {
		String sequenceName = "sequence";
		String first = "first";
		String last = "last";
		String sample_count = "methylated_sample_count";
		String diff = "methylation_percentage";
		String pVal = "p-value";
		String correctedpVal = "BH_corrected_alpha";
		writer.println(sequenceName + "\t" + first + "\t" + last + "\t" + 
				sample_count + "\t"+ diff + "\t" + pVal + "\t" + correctedpVal);
	}
	
	public void printMrs() {
		List<MethylatedRegion> mrsList = methylatedRegions.asList();
		for (MethylatedRegion mr : mrsList) {
			String mrKey = mr.decodeMRToString();
			int calls = mrSampleCounts.get(mrKey);
			mr.setCalls(calls);
			String mrRecord = mr.toString();
			writer.println(mrRecord);
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String refFile = args[0];
		String outPrefix = args[1];
		int swLength = Integer.parseInt(args[2]);
		double initialAlpha = Double.parseDouble(args[3]);
		MethylatedRegionFinder instance = new MethylatedRegionFinder(refFile, swLength, initialAlpha);
		List<String> files = new ArrayList<>();
		for(int f = 4; f < args.length; f++) {
			files.add(args[f]);
		}
		instance.run(files, outPrefix);
	}

}
