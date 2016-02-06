package org.geatools.operation;
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
import  java.io.*;  
import java.util.ArrayList; 
import java.util.List;


public  class  FileOperate  {
  
   public  FileOperate()  {  
   } 
   
   static public List<ArrayList <String>> getMatrixFromFile(String inFile){
  
    List<ArrayList <String>> outList=new ArrayList<ArrayList <String>>();
    ArrayList <String> lineList=new ArrayList <String> ();
	try{    
        BufferedReader br;              
        br = new BufferedReader(new FileReader(inFile));
	    String line;
		String [] itemSplited;
		while(true){
		   line = br.readLine();
	       if (line == null) break;
	       itemSplited=line.split("\t",-1);
		   lineList=new ArrayList <String> ();
		   for(int i=0;i<itemSplited.length;i++){
		      lineList.add(itemSplited[i].trim());
		   }
		   outList.add(lineList);		
		}
		br.close();
		lineList=null;
	}catch(IOException e){
        System.out.println(e);
    }  

    return outList;	
  
   }  
      
   /**  
    *  Save ArrayList<ArrayList<String>> as file  
    *  @param ArrayList<ArrayList<String>> str, outFile
   */  
  
  static public void saveMatrixList(List<ArrayList<String>> str,File outFile){
 
   try{
           String text="";
           BufferedWriter writer=null;
           writer=new BufferedWriter(new FileWriter(outFile));

           for(int i=0;i<str.size();i++){
			  text="";
             for(int j=0;j<str.get(i).size();j++){
               
				if(j<str.get(i).size()-1)
				 text=text+str.get(i).get(j)+"\t";
				else if(j==str.get(i).size()-1)
				 text=text+str.get(i).get(j);
             }
             writer.write(text);
             writer.newLine();
             writer.flush();     
			  
           } 
			
			writer.close();
        
           text=null;
   }catch(IOException e){
           System.out.println(e);
   }     
  }
   
   /**  
     *  Save ArrayList<ArrayList<String>> as file  
     *  @param ArrayList<ArrayList<String>> str, outFile
    */  
   
   static public void saveMatrixList(List<ArrayList<String>> str,String outFile){
  
    try{
            String text="";
            BufferedWriter writer=null;
            writer=new BufferedWriter(new FileWriter(outFile));

            for(int i=0;i<str.size();i++){
			  text="";
              for(int j=0;j<str.get(i).size();j++){
                
				if(j<str.get(i).size()-1)
				 text=text+str.get(i).get(j)+"\t";
				else if(j==str.get(i).size()-1)
				 text=text+str.get(i).get(j);
              }
              writer.write(text);
              writer.newLine();
              writer.flush();     
			  
            } 
			
			writer.close();
         
            text=null;
    }catch(IOException e){
            System.out.println(e);
    }     
   }
   
   static public void saveList(List<String> str, String header, File outFile){
	   
	    try{
	            String text=header;
	            BufferedWriter writer=null;
	            writer=new BufferedWriter(new FileWriter(outFile));
	            if(header!=null && header!=""){
		          writer.write(text);
		          writer.newLine();
		          writer.flush(); 
	            }

	            for(int i=0;i<str.size();i++){
	              text=str.get(i);
	              
	              writer.write(text);
	              writer.newLine();
	              writer.flush();    
	            } 
				
				writer.close();
	            text=null;
	    }catch(IOException e){
	            System.out.println(e);
	    }     
   }
   
   static public void saveList(List<String> str,String header,String outFile){
	   
	    try{
	            String text=header;
	            BufferedWriter writer=null;
	            writer=new BufferedWriter(new FileWriter(outFile));
                if(header!=null && header!=""){
  	              writer.write(text);
  	              writer.newLine();
  	              writer.flush(); 
                }
                
	            for(int i=0;i<str.size();i++){
	              text=str.get(i);
	              
	              writer.write(text);
	              writer.newLine();
	              writer.flush();    				  
	            } 
				
				writer.close();
	            text=null;
	    }catch(IOException e){
	            System.out.println(e);
	    }     
   }
   
   static public List<String> getRowListFromFile(String inFile){
		 
	    List <String> lineList=new ArrayList <String> ();
		try{    
	        BufferedReader br;              
	        br = new BufferedReader(new FileReader(inFile));
		    String line;
		
			while(true){
			   line = br.readLine();
		       if (line == null || line.equals("")) break;

			   lineList.add(line.trim());		
			}
			br.close();
			
		}catch(IOException e){
	        System.out.println(e);
	    }  

	    return lineList;	
	  
   }
  
   static public List<String> combineRowListFromFiles(String inFilesFile, String outFile){
		 
	    List <String> lineList=new ArrayList <String> ();
		try{    
	      BufferedReader br;  
	   	  List<String> seqFileList=FileOperate.getRowListFromFile(inFilesFile);
	      if(outFile==null) outFile=inFilesFile+".rowCombine";
	   	  BufferedWriter writer=null;
          writer=new BufferedWriter(new FileWriter(outFile));
   	      for(String seq:seqFileList){    		
	        br = new BufferedReader(new FileReader(seq));
		    String line;		
			while(true){
			   line = br.readLine();
		       if (line == null || line.equals("")) break;
		       writer.write(line);
		       writer.newLine();
		       writer.flush(); 	         
			}
			br.close();
   	      }
		  writer.close();	
		}catch(IOException e){
	        System.out.println(e);
	    }  

	    return lineList;	
	  
   }
  
   static public List<String> combineRowListFromFiles(List<String> seqFileList,
		  String outFile){
		 
	    List <String> lineList=new ArrayList <String> ();
		try{    
	      BufferedReader br;  
	 
	      if(outFile==null) outFile="rowNamesCombine.txt";
	   	  BufferedWriter writer=null;
         writer=new BufferedWriter(new FileWriter(outFile));
  	      for(String seq:seqFileList){    		
	        br = new BufferedReader(new FileReader(seq));
		    String line;		
			while(true){
			   line = br.readLine();
		       if (line == null || line.equals("")) break;
		       writer.write(line);
		       writer.newLine();
		       writer.flush(); 	         
			}
			br.close();
  	      }
		  writer.close();	
		}catch(IOException e){
	        System.out.println(e);
	    }  

	    return lineList;	
	  
   }

   
   public static String getFileFormat(String file){
		 
	     if(file==null) return null;
		 
	     String format=file.substring(file.lastIndexOf(".")+1,file.length());
		 
		 return format;
   }

   /**  
     *  
     *  @param  folderPath  String
     *  @return  boolean  
     */  
   static public void  newFolder(String  folderPath)  {  
       try  {  
           String  filePath  =  folderPath;  
           filePath  =  filePath.toString();  
           java.io.File  myFilePath  =  new  java.io.File(filePath);  
           if  (!myFilePath.exists())  {  
               myFilePath.mkdir();  
           }  
       }  
       catch  (Exception  e)  {  
           System.out.println("Fail to create folder");  
           e.printStackTrace();  
       }  
   }  

   /**  
     *  
     *  @param  filePathAndName  String 
     *  @param  fileContent  String  
     *  @return  boolean  
     */  
   static public void  newFile(String  filePathAndName,  String  fileContent)  {  

       try  {  
           String  filePath  =  filePathAndName;  
           filePath  =  filePath.toString();  
           File  myFilePath  =  new  File(filePath);  
           if  (!myFilePath.exists())  {  
               myFilePath.createNewFile();  
           }  
           FileWriter  resultFile  =  new  FileWriter(myFilePath);  
           PrintWriter  myFile  =  new  PrintWriter(resultFile);  
           String  strContent  =  fileContent;  
           myFile.println(strContent);  
           resultFile.close();  

       }  
       catch  (Exception  e)  {  
           System.out.println("Fail to create file");  
           e.printStackTrace();  

       }  

   }  

   /**  
     *    
     *  @param  filePathAndName  String 
     *  @param  fileContent  String  
     *  @return  boolean  
   */  
   static public void  delFile(String  filePathAndName)  {  
       try  {  
           String  filePath  =  filePathAndName;  
           filePath  =  filePath.toString();  
           java.io.File  myDelFile  =  new  java.io.File(filePath);  
           myDelFile.delete();  

       }  
       catch  (Exception  e)  {  
           System.out.println("Fail to delete file");  
           e.printStackTrace();  

       }  

   }  

   /**  
     * 
     *  @param  filePathAndName  String  
     *  @param  fileContent  String  
     *  @return  boolean  
   */  
   static public void  delFolder(String  folderPath)  {  
       try  {  
           delAllFile(folderPath);  //
           String  filePath  =  folderPath;  
           filePath  =  filePath.toString();  
           java.io.File  myFilePath  =  new  java.io.File(filePath);  
           myFilePath.delete();  //

       }  
       catch  (Exception  e)  {  
           System.out.println("Fail to delete folder");  
           e.printStackTrace();  

       }  

   }  

   /**  
     *  
     *  @param  path  String  
     */  
   static public void  delAllFile(String  path)  {  
       File  file  =  new  File(path);  
       if  (!file.exists())  {  
           return;  
       }  
       if  (!file.isDirectory())  {  
           return;  
       }  
       String[]  tempList  =  file.list();  
       File  temp  =  null;  
       for  (int  i  =  0;  i  <  tempList.length;  i++)  {  
           if  (path.endsWith(File.separator))  {  
               temp  =  new  File(path  +  tempList[i]);  
           }  
           else  {  
               temp  =  new  File(path  +  File.separator  +  tempList[i]);  
           }  
           if  (temp.isFile())  {  
               temp.delete();  
           }  
           if  (temp.isDirectory())  {  
               delAllFile(path+"/"+  tempList[i]);//
               delFolder(path+"/"+  tempList[i]);//
           }  
       }  
   }  

   /**  
     *  
     *  @param  oldPath  String 
     *  @param  newPath  String 
     *  @return  boolean  
     */  
   static public void  copyFile(String  oldPath,  String  newPath)  {  
       try  {  
           //int  bytesum  =  0;  
           int  byteread  =  0;  
           File  oldfile  =  new  File(oldPath);  
           if  (oldfile.exists())  {  //
               InputStream  inStream  =  new  FileInputStream(oldPath);  //
               FileOutputStream  fs  =  new  FileOutputStream(newPath);  
               byte[]  buffer  =  new  byte[1444];  
               while  ((byteread  =  inStream.read(buffer))  !=  -1)  {  
                   //bytesum  +=  byteread;  //
                   //System.out.println(bytesum);  
                   fs.write(buffer,  0,  byteread);  
               }  
               inStream.close(); 
               fs.close();
           }           
           
       }  
       catch  (Exception  e)  {  
           System.out.println("Fail to copy file");  
           e.printStackTrace();  

       }  

   }  

   /**  
     *  
     *  @param  oldPath  String 
     *  @param  newPath  String    
     *  @return  boolean  
     */  
   static public void  copyFolder(String  oldPath,  String  newPath)  {  

       try  {  
           (new  File(newPath)).mkdirs();  //
           File  a=new  File(oldPath);  
           String[]  file=a.list();  
           File  temp=null;  
           for  (int  i  =  0;  i  <  file.length;  i++)  {  
               if(oldPath.endsWith(File.separator)){  
                   temp=new  File(oldPath+file[i]);  
               }  
               else{  
                   temp=new  File(oldPath+File.separator+file[i]);  
               }  

               if(temp.isFile()){  
                   FileInputStream  input  =  new  FileInputStream(temp);  
                   FileOutputStream  output  =  new  FileOutputStream(newPath  +  "/"  + 
                            (temp.getName()).toString());  
                   byte[]  b  =  new  byte[1024  *  5];  
                   int  len;  
                   while  (  (len  =  input.read(b))  !=  -1)  {  
                       output.write(b,  0,  len);  
                   }  
                   output.flush();  
                   output.close();  
                   input.close();  
               }  
               if(temp.isDirectory()){//
                   copyFolder(oldPath+"/"+file[i],newPath+"/"+file[i]);  
               }  
           }  
       }
       catch  (Exception  e)  {  
           System.out.println("Fail to copy folder");  
           e.printStackTrace();  

       }  

   }  

   /**  
     * 
     *  @param  oldPath  String  
     *  @param  newPath  String  
     */  
   static public void  moveFile(String  oldPath,  String  newPath)  {  
       copyFile(oldPath,  newPath);  
       delFile(oldPath);  

   }  

   /**  
     * 
     *  @param  oldPath  String  
     *  @param  newPath  String  
     */  
   static public void  moveFolder(String  oldPath,  String  newPath)  {  
       copyFolder(oldPath,  newPath);  
       delFolder(oldPath);  

   }

   static public void combineFiles(ArrayList<String> fileList, String outFile){
    if(fileList.size()>0){
	 try{
      String file1=fileList.get(0);
	  String file2="";
	  copyFile(file1, outFile); 
	  
      String line="";
	  BufferedWriter writer=null;
	  writer=new BufferedWriter(new FileWriter(outFile,true));
	  BufferedReader br;  
	  for(int i=1;i<fileList.size();i++){	
   	    file2= fileList.get(i);         
        br = new BufferedReader(new FileReader(file2));
	    while(true){
           line = br.readLine();
	       if (line == null) break;
		   writer.write(line);
		   writer.newLine();
		   writer.flush(); 
		}           		   
		br.close();		
	  
	  }
	  writer.close();
	 }catch(IOException e){
        System.out.println(e);
     }
	}
   }   
} 
