package PostDMRProcessing;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import MethylationDataClasses.MethylatedRegion;

public class JointSamplesDMRAnalyses {
	
	public final static String DEFAULT_CONTROL_SAMPLE_ID = "FC";
	
	private String [] sampleIds;
	private Map<String, Map<String, MethylatedRegion>> dmrs;
	
	public JointSamplesDMRAnalyses(int length) throws IOException {
		sampleIds = new String [length+1];
		dmrs = new LinkedHashMap<>();
	}
	
	public static void main(String[] args) throws IOException  {
		// TODO Auto-generated method stub
		JointSamplesDMRAnalyses instance = new JointSamplesDMRAnalyses(args.length);
		String [] sampleIds = new String [args.length + 1];
		sampleIds[0] = DEFAULT_CONTROL_SAMPLE_ID;
		for(int i = 0; i < args.length; i++){
			String file = args[i];
			boolean firstSample = i == 0;
			instance.readSample(file, firstSample);
			sampleIds[i+1] = file.split("\\.")[0];
		}
		instance.setSampleIds(sampleIds);
		instance.processDMRsToPercentageMatrix();
		instance.printJointRegions();
	}
	private void printJointRegions() throws IOException {
		try(PrintWriter writer = new PrintWriter("jointDmrList.csv")){
			for(Map.Entry<String, Map<String, MethylatedRegion>> entry : dmrs.entrySet()) {
				String region = entry.getKey();
				writer.println(region);
			}
		}
	}
	private void processDMRsToPercentageMatrix() throws IOException {
		try(PrintWriter writer = new PrintWriter("samplePercentageMatrix.tsv")){
			String samplePrefix = "obs.";
			String header = "\t";
			int n = dmrs.keySet().size();
			for(int i = 0; i < n; i++) {
				int nS = i+1;
				if(i < n - 1) header += samplePrefix + nS + "\t";
				else header += samplePrefix + nS;
			}
			writer.println(header);
			for(String sampleId : sampleIds) {
				writer.print(sampleId + "\t");
				for(Map.Entry<String, Map<String, MethylatedRegion>> entry : dmrs.entrySet()) {
					MethylatedRegion current = entry.getValue().get(sampleId);
					if(current == null) writer.print(0 + "\t");
					else writer.print(current.getMethylatedTreatmentPercentage() + "\t");
				}
				writer.println();
			}
		}
	}
	private void readSample(String file, boolean firstSample) throws IOException {
		// TODO Auto-generated method stub
		try(BufferedReader reader = new BufferedReader(new FileReader(file))){
			reader.readLine();
			String sampleId = file.split("\\.")[0];
			String line = reader.readLine();
			while(line != null) {
				String [] elements = line.split("\t");
				String seqName = elements[0];
				int first = Integer.parseInt(elements[1]);
				int last = Integer.parseInt(elements[2]);
				double controlPercentage = Double.parseDouble(elements[3]);
				double treatmentPercentage = Double.parseDouble(elements[4]);
				double pVal = Double.parseDouble(elements[6]);
				String key = encodeKey(seqName, first, last);
				MethylatedRegion treatDmr = new MethylatedRegion(seqName, first, last, controlPercentage, 
						treatmentPercentage, pVal, true, null);
				treatDmr.setSampleId(seqName);
				dmrs.computeIfAbsent(key, k -> new LinkedHashMap<>()).put(sampleId, treatDmr);
				if(firstSample) {
					MethylatedRegion controlDmr = new MethylatedRegion(seqName, first, last, 0, 
							controlPercentage, 0, true, null);
					controlDmr.setSampleId(DEFAULT_CONTROL_SAMPLE_ID);
					dmrs.computeIfAbsent(key, k -> new LinkedHashMap<>()).put(DEFAULT_CONTROL_SAMPLE_ID, controlDmr);
				}
				line = reader.readLine();
			}
		}
	}
	private String encodeKey(String seq, int first, int last) {
		return seq + "," + first  + "," + last;
	}
	public void setSampleIds(String[] sampleIds) {
		this.sampleIds = sampleIds;
	}
}
