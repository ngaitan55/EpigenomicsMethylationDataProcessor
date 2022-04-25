package PostDMRProcessing;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import MethylationDataClasses.MethylatedRegion;
import MethylationDataClasses.MethylationRecord;
import MethylationDataClasses.MethylationSampleFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequenceList;

public class AllReplicatesDMRAnalyses {
	
	public static final double METHYLATION_POS_TRESHOLD = 0.3;
	
	private String [] sampleIds;
	private ReferenceGenome refGenome;
	private GenomicRegionSortedCollection<GenomicRegion> jointDMRRegions;
	private Map<String, Map<String, MethylatedRegion>> dmrs;
	
	public AllReplicatesDMRAnalyses(int nSamples, String refFile) throws IOException {
		this.sampleIds = sampleIds;
		this.refGenome = new ReferenceGenome(refFile);
		this.jointDMRRegions = new GenomicRegionSortedCollection<>(refGenome.getSequencesList());
		dmrs = new LinkedHashMap<>();
	}
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String refFile = args[0];
		AllReplicatesDMRAnalyses instance = new AllReplicatesDMRAnalyses(args.length-2, refFile); 
		String regionsListFile = args[1];
		instance.readDMRRegions(regionsListFile);
		String [] sampleIds = new String[args.length-2];
		for(int i = 2; i < args.length; i++) {
			String file = args[i];
			sampleIds[i - 2] = file.split("\\.")[0];
			instance.readAndComputeRegionsInSample(file, sampleIds[i - 2]);
		}
		instance.setSampleIds(sampleIds);
		instance.processDMRsToPercentageMatrix();
	}
	private void processDMRsToPercentageMatrix() throws IOException {
		// TODO Auto-generated method stub
		try(PrintWriter writer = new PrintWriter("replicatePercentageMatrix.tsv")){
			//String samplePrefix = "obs.";
			String samplePrefix = "";
			String header = "\t";
			int n = dmrs.keySet().size();
			int i = 0;
			/**
			for(int i = 0; i < n; i++) {
				int nS = i+1;
				if(i < n - 1) header += samplePrefix + nS + "\t";
				else header += samplePrefix + nS;
			}**/
			for(String key : dmrs.keySet()) {
				samplePrefix = key;
				if(i < n - 1) header += samplePrefix  + "\t";
				else header += samplePrefix;
				i++;
			}
			writer.println(header);
			for(String sampleId : sampleIds) {
				writer.print(sampleId + "\t");
				for(Map.Entry<String, Map<String, MethylatedRegion>> entry : dmrs.entrySet()) {
					MethylatedRegion current = entry.getValue().get(sampleId);
					writer.print(current.getMethylatedTreatmentPercentage() + "\t");
				}
				writer.println();
			}
		}
	}
	private void readAndComputeRegionsInSample(String file, String sampleId) throws IOException {
		// TODO Auto-generated method stub
		try(MethylationSampleFileReader reader = new MethylationSampleFileReader(file)){
			Iterator<MethylationRecord> it = reader.iterator();
			List<GenomicRegion> regionsList = jointDMRRegions.asList();
			int n = regionsList.size();
			boolean[] visited = new boolean [n];
			int i = 0;
			MethylationRecord mR = it.next();
			int methCount = 0;
			while(i < n) {
				GenomicRegion region = regionsList.get(i);
				int regionLength = region.length();
				int cmp = compare(region, mR);
				//System.out.println(region.getSequenceName() + " " + region.getFirst() + " " + region.getLast());
				//System.out.println(mR.getSequenceName() + " " + mR.getFirst());
				//System.out.println("cmp="+ cmp);
				if(cmp == -1) {
					if(visited[i]) {
						double percentage = ((double) methCount / regionLength)*100;
						computeDMREntry(region, percentage, sampleId);
						methCount = 0;
					}
					else {
						double percentage = 0;
						computeDMREntry(region, percentage, sampleId);
					}
					i++;
				}
				else if(cmp == 0) {
					int methylated = mR.getMethylatedBaseCalls();
					int total = mR.getTotal();
					boolean ans = baseIsMethylated(methylated, total);
					if(ans) methCount++;
					visited[i] = true;
					if(it.hasNext()) mR = it.next();
					else {
						double percentage = ((double) methCount / regionLength)*100;
						computeDMREntry(region, percentage, sampleId);
						methCount = 0;
						i++;
					}
				}
				else {
					if(it.hasNext()) mR = it.next();
					else {
						double percentage = 0;
						computeDMREntry(region, percentage, sampleId);
						i++;
					}
				}
			}
		}
	}
	private void computeDMREntry(GenomicRegion region, double percentage, String sampleId) {
		String seqName = region.getSequenceName();
		int first = region.getFirst();
		int last = region.getLast();
		MethylatedRegion treatDmr = new MethylatedRegion(seqName,
				first, last, 0, percentage, 0, true, null);
		treatDmr.setSampleId(sampleId);
		String key = encodeKey(seqName, first, last);
		dmrs.computeIfAbsent(key, k -> new LinkedHashMap<>()).put(sampleId, treatDmr);
	}
	private void readDMRRegions(String regionsListFile) throws IOException {
		// TODO Auto-generated method stub
		try(BufferedReader reader = new BufferedReader(new FileReader(regionsListFile))){
			String line = reader.readLine();
			while(line != null) {
				String [] elements = line.split(",");
				String seqName = elements[0];
				int first = Integer.parseInt(elements[1]);
				int last = Integer.parseInt(elements[2]);
				jointDMRRegions.add(new GenomicRegionImpl(seqName, first, last));
				line = reader.readLine();
			}
			jointDMRRegions.forceSort();
		}
	}
	public int compare(GenomicRegion sw, MethylationRecord m) {
		QualifiedSequenceList sequences = refGenome.getSequencesList();
		int p1 = sequences .indexOf(sw.getSequenceName());
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
	private String encodeKey(String seq, int first, int last) {
		return seq + "," + first  + "," + last;
	}
	private void setSampleIds(String[] samples) {
		this.sampleIds = samples;
	}
}

