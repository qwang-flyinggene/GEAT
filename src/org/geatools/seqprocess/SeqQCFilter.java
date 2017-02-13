package org.geatools.seqprocess;
/*******************************************************************************
 *  ========================================================================
 *  GEATools : a free Genomic Event Analysis Tools for the Java(tm) platform
 *  ========================================================================
 *
 *  (C) Copyright 2016, by Qi Wang and Contributors.
 *
 *  This file is part of GEATools.
 *
 *  GEATools is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 *  [Oracle and Java are registered trademarks of Oracle and/or its affiliates. 
 *  Other names may be trademarks of their respective owners.]
 *
 *  (C) Copyright 2000-2014, by Original Author and Contributors.
 *
 *  Original Author: Qi Wang;
 *  Contributor(s): 
 *  Changes (from 1-Jan-2016)
 *  ---------------------------------------------------------------------------
 *******************************************************************************/
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.geatools.GEAT;
import org.geatools.data.structure.SeqQual;

public class SeqQCFilter {
	String prinseqDir="utilities/prinseq";
	String outDir;
	String resFasta;
	String resFasta2;
	String resFastq;
	String resFastq2;
	String resQual;
	String resQual2;
	
	boolean removeExactDup=false;	
	int minQualMean=-1;
	int minLen=-1;
	int outFormat=-1;
	
	public boolean isOK=true;
	
	public SeqQCFilter(){
		
	}
	public void setPrinseqDir(String dir){
		prinseqDir=dir;
	}
	
	public String getPrinseqDir(){
		return prinseqDir;
	}
	
	public void setOutDir(String dir){
		outDir=dir;
	}
	
	public String getResFasta(){
		return resFasta;
	}	
	
	public String getResFastq(){
		return resFastq;
	}
	
	public String getResQual(){
		return resQual;
	}
	
	/*
	public String getResFasta2(){
		return resFasta2;
	}
	
	public String getResFastq2(){
		return resFastq2;
	}
	
	public String getResQual2(){
		return resQual2;
	}
	*/
	public void setOpts(String[]opts){
	   if(opts!=null){ 
		 Map<String, List<String>> params=GEAT.getCommandLineParams(opts);

		 //####### set seq minmum length #######
		 if(params.get("min_len")!=null){		  	
		   if(params.get("min_len").size()>0){
			  String tmp=params.get("min_len").get(0).trim();
			  try{
				 minLen=Integer.parseInt(tmp);
				 if(minLen<0){ 
					 System.err.println("Illegal '-min_len' parameter usage :(");
					 return;
				 }	
			  }catch(Exception e){
				 System.err.println("Illegal '-min_len' parameter usage :(");
				 return;
			  } 
		   }
		 }
		 
		 //####### set seq minmum quality mean #######
		 if(params.get("min_qual_mean")!=null){		  	
		   if(params.get("min_qual_mean").size()>0){
			  String tmp=params.get("min_qual_mean").get(0).trim();
			  try{
				 minQualMean=Integer.parseInt(tmp);
				 if(minQualMean<0){ 
					 System.err.println("Illegal '-min_qual_mean' parameter usage :(");
					 return;
				 }	
			  }catch(Exception e){
				 System.err.println("Illegal '-min_qual_mean' parameter usage :(");
				 return;
			  } 
			
		   }
		 }
		 
		 //####### set removal of exact duplication #######
		 if(params.get("de_exact_dup")!=null){	
		   removeExactDup=true;
		   try{
		     if(params.get("de_exact_dup").size()>0){
			     String tmp=params.get("de_exact_dup").get(0).trim();
			 				 
				 if(tmp.equalsIgnoreCase("no") 
						 || tmp.equalsIgnoreCase("N")
						 || tmp.equalsIgnoreCase("false")
						 || tmp.equalsIgnoreCase("F")){ 

					 removeExactDup=false;
				 }				 
		      }
		   }catch(Exception e){
			  System.err.println("Illegal '-de_exact_dup' parameter usage :(");
			  return;
		   } 
		 }
		 
		//####### set out_format #######
		 if(params.get("out_format")!=null){		  	
		   if(params.get("out_format").size()>0){
			  String tmp=params.get("out_format").get(0).trim();			
			  try{
				 outFormat=Integer.parseInt(tmp);
				 if(outFormat<=0){ 
					 System.err.println("Illegal '-out_format' parameter usage :(");
					 return;
				 }	
			  }catch(Exception e){
				 System.err.println("Illegal '-out_format' parameter usage :(");
				 return;
			  } 
			  
		   }
		 }
	   }
	}
	
	public List<SeqQual> getSeqQual(String seqQualFile){
		
		if(!new File(seqQualFile).exists()) return null;
		List<SeqQual> seqQualList=new ArrayList<SeqQual>();
		SeqQual seqQual;
		int [] seqQualArray;
		int seqNum=0;
		try{    
	        BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqQualFile));
		    String line;
			String seqNameLine;
			String seqQualLine;
			String seqName;
			String [] itemSplited;				
			seqNum=0;
			line = br.readLine();				
			while(true){           
		       if (line == null) break;
		       line=line.trim();
			   if(line.indexOf(">")==0){
			     seqNum=seqNum+1;
			     seqNameLine=line;				 
				 itemSplited=seqNameLine.split("\\s+");
				 seqName=itemSplited[0].trim();
				 seqName=seqName.substring(seqName.indexOf(">")+1,seqName.length());
				 
				 line=br.readLine();
				 seqQualLine="";
				 if (line == null) break;
				 while(line.indexOf(">")<0){
				   seqQualLine = seqQualLine+line.trim()+" ";
				   line=br.readLine();
				   if (line == null) break;
				 }
				 seqQualLine=seqQualLine.trim();
				 itemSplited=seqQualLine.split("\\s+");
				 seqQualArray=new int[itemSplited.length];
				 for(int i=0;i<seqQualArray.length;i++){				
					seqQualArray[i]=Integer.parseInt(itemSplited[i].trim());
				 }
				 itemSplited=null;
				 seqQual=new SeqQual();
				 seqQual.seqName=seqName;
				 seqQual.seqQual=seqQualArray;
				 seqQualList.add(seqQual);
	             seqQual=null;		             
			   }				 
			}           		   
			br.close();
		
		}catch(IOException e){
	        System.out.println(e);
	    }

		return seqQualList;
	}
	
	public void prinseqQC(String fastqFile){
		   
		   isOK=true;
		   
		   File f=new File(fastqFile);
		   if(!f.exists()){
			  isOK=false;
			  f=null; 
			  return;
		   }
		   
		   String seqFile=fastqFile;
		   if(outDir==null) outDir=seqFile.substring(0,seqFile.lastIndexOf("/"));		   
		   		 
		   String seqName=seqFile.substring(
				            seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
				          );
		   String seqFullName=outDir+"/"+seqName;
		   	   
		 
		   Process process;
		   try {
			 List<String> cmdList=new ArrayList<String>();
			 List<String> cmdStrList=new ArrayList<String>();		   	  
		   	 
		   	 if(minQualMean>0 || minLen>0){
		   		 
		   		 String cmd="perl "+prinseqDir+"/prinseq-lite.pl -fastq "+seqFile;
			   	 
			   	 if(minQualMean>0) cmd=cmd+" -min_qual_mean "+minQualMean;
			   	 
			   	 if(minLen>0) cmd=cmd+" -min_len "+minLen;
			   	 
			   	 if(removeExactDup) cmd=cmd+" -out_format 3";
			   	 else if(outFormat>0) cmd=cmd+" -out_format "+outFormat;			   	 
			
			   	 cmd=cmd+" -out_good "+seqFullName+".goodQual -out_bad "+seqFullName+".badQual";	
					   	 
			   	 cmdList.add(cmd);
			     cmdStrList.add("Filtering sequence with bad quality......"); 
			   	 
			     if(outFormat<=2 && !removeExactDup) seqFile=seqFullName+".goodQual.fasta";
			     else seqFile=seqFullName+".goodQual.fastq";
			     
			     //seqName=seqFile.substring(0,seqFile.lastIndexOf("."));
			     seqName=seqFile.substring(
				            seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
				          );
			     seqFullName=outDir+"/"+seqName;
		   	 }
		   	 
		   	 if(removeExactDup){
		   		
		   		 String cmd="perl "+prinseqDir+"/prinseq-lite.pl -fastq "+seqFile
		   	    		    +" -derep 14 -derep_min 2 ";
			     if(outFormat>0) cmd=cmd+"-out_format "+outFormat;
			     cmd=cmd+" -out_good "+seqFullName+".goodDup -out_bad "+seqFullName+".badDup";
			  
			     cmdList.add(cmd); 
			     cmdStrList.add("Filtering sequence with exact duplication......");
			     
			     if(outFormat<=2) seqFile=seqFullName+".goodDup.fasta";
			     else if(outFormat<=5) seqFile=seqFullName+".goodDup.fastq";
			     
			     //seqName=seqFile.substring(0,seqFile.lastIndexOf("."));
			     seqName=seqFile.substring(
				            seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
				          );
			     seqFullName=outDir+"/"+seqName;
		   	 }
		   	 
		  	 if(outFormat>0 && minQualMean<0 && minLen<0 && !removeExactDup){
		  	
		  		 String cmd="perl "+prinseqDir+"/prinseq-lite.pl -fastq "+seqFile;
			   	 
			   	 cmd=cmd+" -out_format "+outFormat;			   	 
			
			   	 cmd=cmd+" -out_good "+seqFullName+".Out";	
					   	 
			   	 cmdList.add(cmd);
			     cmdStrList.add("Outputing Sequences......"); 
			   	 
			     if(outFormat<=2) seqFile=seqFullName+".Out.fasta";
			     else if(outFormat<=5) seqFile=seqFullName+".Out.fastq";
			     
			     //seqName=seqFile.substring(0,seqFile.lastIndexOf("."));
			     seqName=seqFile.substring(
				            seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
				          );
			     seqFullName=outDir+"/"+seqName;
		   	 }
	
			
			 for(int i=0;i<cmdList.size();i++){
			   
			   System.out.println("Step "+(i+1)+": "+cmdStrList.get(i));
			   
			   process = Runtime.getRuntime().exec(cmdList.get(i));
			   process.waitFor();
			   if(process.exitValue() == 0) {
			     System.out.println("Success!");
			   }else {
				 isOK=false;
			     System.err.println("Failure!!!");
			     return;
			   }
			   
			 }
			 
			 if(isOK){
				if(outFormat==1 || outFormat==2 || outFormat==4 || outFormat==5){
			       resFasta=seqFullName+".fasta";
				}
				
				if(outFormat==3 || outFormat==4 || outFormat==5){
				   resFastq=seqFullName+".fastq";
				}
				
				if(outFormat==2 || outFormat==5){
				   resQual=seqFullName+".qual";
				}
			 }else{
			    resFasta=null;
			    resFastq=null;
			 }
		   } catch(Exception e) {
			    System.out.println("Exception: "+ e.toString());
		   }
	}
	
}
