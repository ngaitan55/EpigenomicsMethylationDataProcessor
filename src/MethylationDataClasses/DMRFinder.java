package MethylationDataClasses;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.function.Function;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.ReferenceGenome;
import ngsep.math.FisherExactTest;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class DMRFinder {
	public static final double METHYLATION_POS_TRESHOLD = 0.3;
	public static final int SMOOTHING_ALPHA = 1;
	
	private double significanceLevel;
	private int slidingWindowSize;
	private String [] customSampleIds;
	private PrintWriter writer = null;
	private ReferenceGenome refGenome;
	private QualifiedSequenceList sequences;
	private GenomicRegionSortedCollection<MethylatedRegion> dmrs;
	private int nHypotheses;
	
	public DMRFinder(String genomeFile, int swLength, double initialAlpha) throws IOException {
		nHypotheses = 0;
		significanceLevel = initialAlpha;
		slidingWindowSize = swLength;
		refGenome = new ReferenceGenome(genomeFile);
		sequences = refGenome.getSequencesList();
		dmrs = new GenomicRegionSortedCollection<>(sequences);
	}
	private void run(List<String> files, String outPrefix) throws IOException {
		// TODO Auto-generated method stub
		assignSampleNames(files);
		writer = new PrintWriter(outPrefix + ".DMR.tsv");
		printHeaderFields();
		readAndProcessSamples(files);
		System.out.println("Applying correction");
		applyCorrection();
		System.out.println("Printing output file");
		printDmrs();
		System.out.println("Printing examples file");
		printExamples(outPrefix + ".examples.tsv");
		System.out.println("Finished");
	}
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		List<String> files = new ArrayList<>();
		String file = args[0];
		String outPrefix = args[1];
		int swLength = Integer.parseInt(args[2]);
		double initialAlpha = Double.parseDouble(args[3]);
		DMRFinder instance = new DMRFinder(file, swLength, initialAlpha);
		for (int i = 4; i < args.length; i++) {
			files.add(args[i]);
			}	
		instance.run(files, outPrefix);
		}
		
	public void readAndProcessSamples(List<String> files) throws IOException{
		List<MethylationSampleFileReader> readers = new ArrayList<>();
		List<Iterator<MethylationRecord>> iterators = new ArrayList<>();
		for (String file:files) {
			MethylationSampleFileReader reader = new MethylationSampleFileReader(file);
			readers.add(reader);
			iterators.add(reader.iterator());
		}
		int nSamples = files.size();
		int nReplicates = nSamples/2;
		int currentSWBegin = 0;
		int currentSWLimit = slidingWindowSize;
		int methylatedSumCS = 0;
		int methylatedSumTS = 0;
		int unmethylatedSumCS = 0;
		int unmethylatedSumTS = 0;
		boolean someHasNext = true;
		for(QualifiedSequence seq : sequences) {
			int seqSize = seq.getLength();
			int nWindows = seqSize/slidingWindowSize;
			int lastWindowSize = seqSize%slidingWindowSize;
			String seqName = seq.getName();
			System.out.println("Processing sequence: " + seqName);
			boolean [] iteratorOnHold = new boolean[nSamples];
			MethylationRecord [] currentNext = new MethylationRecord[nSamples];
			for(int i = 0; i < nWindows; i++) {
				currentSWLimit = i + 1 == nWindows ? lastWindowSize : currentSWLimit;
				GenomicRegion slidingWindow = new GenomicRegionImpl(seqName, currentSWBegin, currentSWLimit);
				SortedMap<Integer, MethylationRecord[] > dmrRecords = new TreeMap<>();
				for(int r = 0; r < nSamples; r++){
					Iterator<MethylationRecord> it = iterators.get(r);
					while(it.hasNext()) {
						if(!iteratorOnHold[r]) {
							currentNext[r] = it.next();
						}
						int cmp = compare(slidingWindow, currentNext[r]);
						if(cmp == -1) {
							iteratorOnHold[r] = true;
							break;
						}
						else if(cmp == 1) {
							iteratorOnHold[r] = false;
							if(!it.hasNext()) 
								break;
						}
						else {
							dmrRecords.computeIfAbsent(currentNext[r].getFirst(),
								k -> new MethylationRecord[nSamples])[r] = currentNext[r];
							int methylated = currentNext[r].getMethylatedBaseCalls();
							int total = currentNext[r].getTotal();
							boolean ans = baseIsMethylated(methylated, total);
							if(r < nReplicates) {
								if(ans) methylatedSumCS ++;
								else unmethylatedSumCS ++;
							}else {
								if(ans) methylatedSumTS ++;
								else unmethylatedSumTS ++;
							}
							iteratorOnHold[r] = false;
						}
					}
				}
				testDifferentialMethylation(nReplicates, methylatedSumCS, unmethylatedSumCS, methylatedSumTS, unmethylatedSumTS,
						slidingWindow, dmrRecords);
				methylatedSumCS = 0;
				unmethylatedSumCS = 0;
				methylatedSumTS = 0;
				unmethylatedSumTS = 0;
				currentSWBegin += slidingWindowSize;
				currentSWLimit += slidingWindowSize;
			}
			currentSWBegin = 0;
			currentSWLimit = slidingWindowSize;
			for(Iterator<MethylationRecord> it :  iterators) {
				someHasNext = it.hasNext();
			}
			System.out.println("Processed sequence: " + seqName);
			if(!someHasNext) break;
		}
		for(MethylationSampleFileReader reader : readers) {
			reader.close();
		}
		dmrs.forceSort();
	}
	
	private void testDifferentialMethylation(int nReplicates, int methylatedSumCS, int unmethylatedSumCS, int methylatedSumTS,
			int unmethylatedSumTS, GenomicRegion slidingWindow, SortedMap<Integer, MethylationRecord[] > dmrRecords) {
		if((methylatedSumCS + unmethylatedSumCS + 
				methylatedSumTS + unmethylatedSumTS) == 0) return;
		String seqName = slidingWindow.getSequenceName();
		int first = slidingWindow.getFirst();
		int last = slidingWindow.getLast();
		int length = slidingWindow.length()*nReplicates;
		smoothContingencyMatrix(methylatedSumCS, unmethylatedSumCS, methylatedSumTS, unmethylatedSumTS, SMOOTHING_ALPHA);
		double controlPerc = ((double) methylatedSumCS/(length))*100;
		double treatPerc = ((double) methylatedSumTS/(length))*100;
		//FisherExact ftInstance = new FisherExact(maxSum);
		double pValue = FisherExactTest.calculatePValue(methylatedSumCS, methylatedSumTS, unmethylatedSumCS, unmethylatedSumTS);
		//double pValue =  ftInstance.getP(methylatedSumCS, methylatedSumTS, unmethylatedSumCS, unmethylatedSumTS);
		nHypotheses++;
		if(pValue <= significanceLevel) {
			MethylatedRegion dmr = new MethylatedRegion(seqName, first, last, controlPerc, 
					treatPerc, pValue, true, dmrRecords);
			dmrs.add(dmr);
		}
	}
	
	public int compare(GenomicRegion sw, MethylationRecord m) {
		int p1 = sequences.indexOf(sw.getSequenceName());
		int p2 = sequences.indexOf(m.getSequenceName());
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
	
	public void smoothContingencyMatrix(int methylatedSumCS, int unmethylatedSumCS, int methylatedSumTS,
			int unmethylatedSumTS, int alpha) {
		methylatedSumCS+=alpha;
		unmethylatedSumCS+=alpha;
		methylatedSumTS+=alpha;
		unmethylatedSumTS+=alpha;
	}
	
	public void applyCorrection () {
		//double m = (double) dmrs.size();
		double m = (double) nHypotheses;
		double alpha =  significanceLevel / m;
		List<MethylatedRegion> dmrsList = dmrs.asList();
		List<MethylatedRegion> toRemove = new ArrayList<>();
		for (MethylatedRegion dmr : dmrsList) {
			if(dmr.getpValue() > alpha) toRemove.add(dmr);
			else dmr.setCorrectedPvalue(dmr.getpValue() * m);
		}
		dmrs.removeAll(toRemove);
		dmrs.forceSort();
	}
	
	public void assignSampleNames(List<String> files) {
		customSampleIds = new String[files.size()];
		int nSamples = files.size();
		int nReps = nSamples/2;
		int repId = 1;
		for(int r = 0; r < nSamples; r++) {
			if(r == nReps) repId = 1;
			customSampleIds[r] = r < nReps ? "C" + repId : "T" + repId;
			repId++;
		}
	}
	
	public void printHeaderFields() {
		String sequenceName = "sequence";
		String first = "first";
		String last = "last";
		String controlPerc = "control_meth_percentage";
		String treatPerc = "treatment_meth_percentage";
		String diff = "percentage_difference";
		String pVal = "p-value";
		String correctedpVal = "corrected_p-value";
		writer.println(sequenceName + "\t" + first + "\t" + last + "\t" + 
				controlPerc + "\t" + treatPerc + "\t" + diff + "\t" + pVal + "\t" + correctedpVal);
	}
	
	public void printDmrs() {
		List<MethylatedRegion> dmrsList = dmrs.asList();
		for (MethylatedRegion dmr : dmrsList) {
			String dmrRecord = dmr.toString();
			writer.println(dmrRecord);
		}
		writer.close();
	}
	
	//Provides the most and least different DMRs based on corrected pValue
	private void printExamples(String examplesFileName) throws IOException {
		PrintWriter examplesWriter = new PrintWriter(examplesFileName);
		List<MethylatedRegion> dmrsList = dmrs.asList();
		Collections.sort(dmrsList, Comparator.comparingDouble((v) -> v.getCorrectedPvalue()));
		printDMRAsTable(dmrsList.get(dmrsList.size() - 1), examplesWriter);
		examplesWriter.println("#-------------------------------------------------------------------------------------------#");
		printDMRAsTable(dmrsList.get(0), examplesWriter);
		examplesWriter.close();
	}
	private void printDMRAsTable(MethylatedRegion dmr, PrintWriter examplesWriter) {
		SortedMap<Integer, MethylationRecord[]> methylatedBases = dmr.getMethylatedBases();
		int firstBase = methylatedBases.firstKey();
		int lastBase = methylatedBases.lastKey();
		int nSamples = customSampleIds.length;
		String exampleHeaders = "pos" + "\t" + "meth_percentage" + "\t" + "meth_calls" + "\t" + "total_calls" + "\t" + "sample";
		examplesWriter.println(exampleHeaders);
		for(int i = firstBase; i <= lastBase; i++) {
			if(methylatedBases.containsKey(i)) {
				MethylationRecord[] baseRecords = methylatedBases.get(i);
				for(int j = 0; j < nSamples; j++) {
					if(baseRecords[j] == null) {
						String genericForNotCalledBase = i + "\t" + 0 + "\t" 
								+ 0 + "\t" + 0 + "\t" + customSampleIds[j];
						examplesWriter.println(genericForNotCalledBase);
					}else {
						MethylationRecord baseRecord = baseRecords[j];
						String recordForCalledBase = i + "\t" + (baseRecord.getMethCallPercentage()*100) + "\t" 
								+ baseRecord.getMethylatedBaseCalls() + "\t" + baseRecord.getTotal() + 
								"\t" + customSampleIds[j];
						examplesWriter.println(recordForCalledBase);
					}
				}
			}
			else {
				for(int j = 0; j < nSamples; j++) {
					String genericForNotCalledBase = i + "\t" + 0 + "\t" 
							+ 0 + "\t" + 0 + "\t" + customSampleIds[j];
					examplesWriter.println(genericForNotCalledBase);
				}
			}
		}
	}
}
