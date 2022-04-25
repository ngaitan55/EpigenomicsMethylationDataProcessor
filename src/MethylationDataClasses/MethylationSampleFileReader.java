package MethylationDataClasses;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.SortedMap;

public class MethylationSampleFileReader implements Iterable<MethylationRecord>,Closeable {
	
	private String sampleId;
	private boolean isReplicate;
	private String replicateId;	
	private Iterator<MethylationRecord> currentIterator = null;
	private BufferedReader in;
	
	public MethylationSampleFileReader(String file) throws IOException {
		init(file);
	}
	@Override
	public Iterator<MethylationRecord> iterator() {
		// TODO Auto-generated method stub
		if(in == null) {
			 throw new IllegalStateException("File reader is closed");
		}
		if(currentIterator != null) {
			throw new IllegalStateException("Iteration in progress");
		}
		currentIterator = new MethylationFileIterator();
		return currentIterator;
	}
	
	private MethylationRecord loadMethylationRecord (String line) {
		String [] items = line.split("\t");
		String sequenceName = items[0];
		int first = Integer.parseInt( items[1]);
		int calls = Integer.parseInt(items[3]);
		int total = Integer.parseInt(items[2]);
		return new MethylationRecord(sequenceName, first, calls, total);
	}
	
	public void init (String file) throws IOException {
		InputStream stream = new FileInputStream(file);
		sampleId = file;
		in = new BufferedReader(new InputStreamReader(stream)); 
	}
	public String getSampleId() {
		return sampleId;
	}
	@Override
	public void close() throws IOException {
		// TODO Auto-generated method stub
		in.close();
	}
	
	private class MethylationFileIterator implements Iterator<MethylationRecord>{
		private MethylationRecord nextRecord;
		
		public MethylationFileIterator() {
			nextRecord = loadRecord();
		}
		@Override
		public MethylationRecord next() {
			if(nextRecord==null) throw new NoSuchElementException();
			MethylationRecord answer = nextRecord;
			nextRecord = loadRecord();
			return answer;
		}
		private MethylationRecord loadRecord() {
			String line;
			while(true) {
				try {
					line = in.readLine();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
				if(line==null) return null;
				MethylationRecord answer = loadMethylationRecord(line);
				if(answer !=null) return answer;
			} 
		}
		@Override
		public boolean hasNext() {
			return nextRecord!=null;
		}
	}
	
	
}
