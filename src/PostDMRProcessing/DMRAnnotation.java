package PostDMRProcessing;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.TranscriptSegment;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;

public class DMRAnnotation {
	
	public static final int SPAN_TO_ANNOTATE = 3000;
	
	public static final String UPSTREAM = "upstream";
	public static final String DOWNSTREAM = "downstream";
	public static final String INTERSECTING = "intersecting";

	public static final String GENE = "gene";
	public static final String CDS = "CDS";
	public static final String R_5P_UTR = "five_prime_UTR";
	public static final String R_3P_UTR = "three_prime_UTR";
	public static final String NC_RNA = "ncRNA";
	
	private ReferenceGenome genome;
	private Transcriptome transcriptome;
	private GenomicRegionSortedCollection<GenomicRegion> dmrs;
	
	public DMRAnnotation(String refFile) throws IOException {
		genome = new ReferenceGenome(refFile);
		dmrs = new GenomicRegionSortedCollection<>(genome.getSequencesList());
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String refFile = args[0];
		String transcriptomeFile = args[1];
		String dmrFile = args[2];
		DMRAnnotation inst = new DMRAnnotation(refFile);
		inst.loadTranscriptome(transcriptomeFile);
		inst.loadDMRs(dmrFile);
		inst.annotateRegionsAndPrint();
	}
	
	public void annotateRegionsAndPrint() throws IOException {
		PrintWriter writer = new PrintWriter("mrRegionsAnnotated.bed");
		String header = "dmr_sequence" + "\t" + "dmr_first" + "\t" + "dmr_last" + "\t" + "annot_sequence" + "\t" +
				"annot_first" + "\t" + "annot_last" + "\t" + "annot_type" + "\t" + "relative_position" + "\t" + "annot_id";
		writer.println(header);
		dmrs.forceSort();
		for(GenomicRegion dmr : dmrs) {
			int upstreamLimit = dmr.getFirst() - SPAN_TO_ANNOTATE;
			int downstreamLimit = dmr.getLast() + SPAN_TO_ANNOTATE;
			GenomicRegionSortedCollection<Transcript> dmrTranscripts = transcriptome.getTranscripts(dmr.getSequenceName(),
					upstreamLimit, downstreamLimit);
			List<Transcript> dmrTranscriptsList = dmrTranscripts.asList();
			for(Transcript transcript : dmrTranscriptsList) {
				printRegionAnnotation(dmr, transcript, writer); 
			}
		}
		writer.close();
	}
	
	public void printRegionAnnotation(GenomicRegion dmr, Transcript transcript, PrintWriter writer) {
		Gene gene = transcript.getGene();
		List<TranscriptSegment> coveredSegments = transcript.getTranscriptSegmentsByAbsolute(dmr.getFirst(), dmr.getLast());
		String position =  determineDmrGeneRelativePosition(dmr, gene);
		String dmrInfo = dmr.getSequenceName() + "\t" + dmr.getFirst() + "\t" + dmr.getLast() + "\t";
		String spanningGeneInfo = gene.getSequenceName() + "\t" + gene.getFirst() + "\t" + gene.getLast() + "\t" + GENE +
				"\t" + position + "\t" + gene.getId();
		String geneLine = dmrInfo + spanningGeneInfo;
		writer.println(geneLine);
		for(TranscriptSegment segment : coveredSegments) {
			String spanningSegmentInfo = segment.getSequenceName() + "\t" + segment.getFirst() + "\t" + segment.getLast() + "\t" +
					translateSegmentType(segment) + "\t" + INTERSECTING + "\t" + segment.getTranscript().getId();
			writer.println(dmrInfo + spanningSegmentInfo);
		}
	}

	private String determineDmrGeneRelativePosition(GenomicRegion dmr, Gene gene) {
		// TODO Auto-generated method stub
		if(gene.getFirst()<=dmr.getFirst() && gene.getLast() >= dmr.getLast() ) {
			return INTERSECTING;
		}
		else if(dmr.getLast() < gene.getFirst()) {
			return gene.isPositiveStrand() ? UPSTREAM : DOWNSTREAM;
		}
		else if(dmr.getFirst() > gene.getLast()) {
			return gene.isPositiveStrand() ? DOWNSTREAM : UPSTREAM;
		}
		return "";
	}

	public String translateSegmentType(TranscriptSegment segment) {
		switch(segment.getStatus()) {
		case 0:
			return CDS;
		case 1:
			return R_5P_UTR;
		case 2:
			return R_3P_UTR;
		case 3:
			return NC_RNA;
		}
		return "";
	}
	
	public void loadTranscriptome(String transcriptomeFile) throws IOException {
		QualifiedSequenceList sequenceNames = genome.getSequencesList();
		GFF3TranscriptomeHandler gff3Handler = new GFF3TranscriptomeHandler(sequenceNames);
		transcriptome = gff3Handler.loadMap(transcriptomeFile); 		
	}
	
	public void loadDMRs(String dmrFile) throws IOException {
		// TODO Auto-generated method stub
		try(BufferedReader reader = new BufferedReader(new FileReader(dmrFile))){
			String line = reader.readLine();
			while(line != null) {
				String [] elements = line.split("\t");
				String seqName = elements[0];
				int first = Integer.parseInt(elements[1]);
				int last = Integer.parseInt(elements[2]);
				dmrs.add(new GenomicRegionImpl(seqName, first, last));
				line = reader.readLine();
			}
		}
	}
}

