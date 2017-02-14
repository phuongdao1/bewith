package BeWithFramework;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class BeWithFramework {
	static HashMap<Integer, String> id2n = new HashMap<Integer, String>();
	static HashMap<String, Integer> n2id = new HashMap<String, Integer>();
	static int MAXMODULESIZE = 10;
	static int numClus = 5;
	static double EDGETHRES = 0.90;

	private static String var(String name, int... args) {
		String ret = name;
		for (int i : args)
			ret = ret + "_" + i;
		return ret;
	}

	private static String coeff(float f) {
		if (f > 0.000001)
			return " + " + f;
		else
			return " - " + (Math.abs(f));
	}

	private static double min(double a, double b) {
		if (a > b)
			return b;
		else
			return a;
	}
	
	public static void howtouse(){
		System.out.println();
		System.out.println("BeWith v0.1: A Between-Within Clustering Framework to Discover Relationships between Cancer Modules");
		System.out.println("             via Integrated Analysis of Mutual Exclusivity, Co-occurrence and Functional Interactions");
		System.out.println("\n");
		
		System.out.println("Usage: java -jar BeWithFramework.jar -m METHOD -k NUMBER_OF_CLUSTERS -o OUTPUT_PREFIX [-NETWORK_TYPE1 INPUT_FILE1] [-NETWORK_TYPE2 INPUT_FILE2]...");
		System.out.println();
		
		System.out.println("METHOD:              bemewithfun, to ensure mutual exclusivity of mutations between different modules");
		System.out.println("                      and functional similarity of genes within modules");
		System.out.println("                     bemewithco, to ensure mutual exclusivity between modules");
		System.out.println("                      and co-occurrence of mutations in genes within modules");
		System.out.println("                     becowithmefun, to ensure co-occurrence between modules while enforcing mutual exclusivity");
		System.out.println("                      and functional interactions within modules");
		System.out.println();
		
		System.out.println("NUMBER OF CLUSTERS:  at least 2 for bemewithfun and bemewithco, currently always set at 2 for becowithfun");
		System.out.println();
		
		System.out.println("OUTPUT_PREFIX:	     the prefix of the result and intermediate files");
		System.out.println();
		
		System.out.println("NETWORK_TYPE:        fun, for a functional or protein interaction network");
		System.out.println("                     me, for a network built from mutually exclusive mutations");
		System.out.println("                     co, for a network built from co-occurring mutations");
		System.out.println();
			
		System.out.println("Example: java -jar BeWithFramework.jar -m bemewithfun -o BRCA_bemewithfun -k 7 -fun data/HumanStringNet.txt -me data/BRCA_me_WESME.txt");
		System.out.println();
		
		//System.out.println("For more reference: ");
		System.out.println();
	}

	public static void main(String[] args) {

		int unselectedCluster = 0;	
		int n = 0;
		int i1, i2;
		double w;
		int i, j, k;
		int u, v;
		int cid;
		String method="";
		String outputPrefix="";
		String funNetworkFile="";
		String meNetworkFile="";
		String coNetworkFile="";
		String outputLPFile="";
				
		if ((args.length%2)!=0){
			System.out.println("ERROR: Invalid passing instructions.");
			howtouse();
			System.exit(1);			
		}
		else if (args.length==0){
			howtouse();
			System.exit(0);	
		}
		
		//	Parsing passing instructions
		for (i=0;i<args.length;i++)
		if (i%2==0){
			String arg=args[i].toLowerCase();
			switch (arg){
				case "-m":
					method=args[i+1].toLowerCase();
					break;
				case "-me":
					meNetworkFile=args[i+1];
					break;
				case "-co":
					coNetworkFile=args[i+1];
					break;
				case "-fun":
					funNetworkFile=args[i+1];
					break;
				case "-o":
					outputPrefix=args[i+1];
					break;
				case "-k":
					numClus = Integer.parseInt(args[i+1])+1;	// actual number of clusters, the cluster of unselected genes is k+1
					break;
				default:
					System.out.println("ERROR: Invalid passing instructions.");
					howtouse();
					System.exit(1);
			}
		}
		
		if (!(method.equals("bemewithfun")||method.equals("bemewithco")||method.equals("becowithmefun"))){
			System.out.println("ERROR: Invalid input method");
			howtouse();
			System.exit(1);	
		}
		
		if (outputPrefix.isEmpty()){
			System.out.println("ERROR: Output name should not be an empty string");
			howtouse();
			System.exit(1);	
		}
		
		if (numClus<=2){
			System.out.println("ERROR: Number of clusters has be at least 2.");
			howtouse();
			System.exit(1);			
		}
		
		if ((funNetworkFile.isEmpty())&&(method.equals("bemewithfun")||method.equals("becowithmefun"))){
			System.out.println("ERROR: Users need to specify a functional or a protein interaction network");
			howtouse();
			System.exit(1);				
		}
		
		if ((coNetworkFile.isEmpty())&&(method.equals("bemewithco")||method.equals("becowithmefun"))){
			System.out.println("ERROR: Users need to specify a network built from co-occurring mutations");
			howtouse();
			System.exit(1);				
		}
		
		if (meNetworkFile.isEmpty()){
			System.out.println("ERROR: Users need to specify a network built from mutually exclusive mutations");
			howtouse();
			System.exit(1);				
		}
		
			
		outputLPFile=outputPrefix+"_"+(numClus-1)+"_modules.lp";
		String line;
		
		HashSet<Integer>[] UG = null;
		float[][] G = null;
		ArrayList<Edge>[] AGL = null;

		float[][] MeG = null;
		ArrayList<Edge>[] AMeGL = null;

		float[][] CoG = null;
		ArrayList<Edge>[] ACoGL = null;

		try {
			// read over the network file to get the number of ids
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(meNetworkFile)));
			//br.readLine();
			while ((line = br.readLine()) != null) {
				line = line.trim();
				String[] tokens = line.split("[\t ]");

				if (!n2id.containsKey(tokens[0])) {
					n2id.put(tokens[0], n);
					id2n.put(n, tokens[0]);
					n++;
				}
				if (!n2id.containsKey(tokens[1])) {
					n2id.put(tokens[1], n);
					id2n.put(n, tokens[1]);
					n++;
				}
			}
			br.close();
			
			if (!coNetworkFile.isEmpty()){
				br = new BufferedReader(new InputStreamReader(new FileInputStream(coNetworkFile)));
				//br.readLine();
				while ((line = br.readLine()) != null) {
					line = line.trim();
					String[] tokens = line.split("[\t ]");
	
					if (!n2id.containsKey(tokens[0])) {
						n2id.put(tokens[0], n);
						id2n.put(n, tokens[0]);
						n++;
					}
					if (!n2id.containsKey(tokens[1])) {
						n2id.put(tokens[1], n);
						id2n.put(n, tokens[1]);
						n++;
					}
				}
				br.close();
			}
			
			
			// preallocate adjacency matrix and list
			G = new float[n][];
			MeG = new float[n][];
			CoG = new float[n][];
			AGL = new ArrayList[n];
			AMeGL = new ArrayList[n];
			ACoGL = new ArrayList[n];
			UG = new HashSet[n];
			for (i = 0; i < n; i++) {
				UG[i] = new HashSet<Integer>();
				G[i] = new float[n];
				MeG[i] = new float[n];
				CoG[i] = new float[n];
				AGL[i] = new ArrayList<Edge>();
				AMeGL[i] = new ArrayList<Edge>();
				ACoGL[i] = new ArrayList<Edge>();
				for (j = 0; j < n; j++) {
					G[i][j] = 0.0f;
					CoG[i][j] = 0.0f;
					MeG[i][j] = 0.0f;
				}
			}

			// read the mutual exclusive graph input file again and fill out
			// adjacent matrix and list
			br = new BufferedReader(new InputStreamReader(new FileInputStream(meNetworkFile)));
			//br.readLine();
			while ((line = br.readLine()) != null) {
				line = line.trim();
				String[] tokens = line.split("[\t ]");

				if (n2id.containsKey(tokens[0]) && n2id.containsKey(tokens[1])) {
					i1 = n2id.get(tokens[0]);
					i2 = n2id.get(tokens[1]);
					if (Float.parseFloat(tokens[2]) <= 0.05) {

						if (Double.parseDouble(tokens[2]) >= 1e-6)
							w = (- Math.log(Double.parseDouble(tokens[2])) / Math.log(100.0));
						else
							w = 3.0;
						

						if (MeG[i1][i2] < 0.000001) {
							UG[i1].add(i2);
							UG[i2].add(i1);
							MeG[i1][i2] = (float) w;
							MeG[i2][i1] = MeG[i1][i2];
							AMeGL[i1].add(new Edge(i2, MeG[i1][i2]));
							AMeGL[i2].add(new Edge(i1, MeG[i1][i2]));
						}
					}
				}
			}
			br.close();

			// read the co-occuring graph input file again and fill out
			// adjacent matrix and list
			if (!coNetworkFile.isEmpty()){
				br = new BufferedReader(new InputStreamReader(new FileInputStream(coNetworkFile)));
				//br.readLine();
				while ((line = br.readLine()) != null) {
					line = line.trim();
					String[] tokens = line.split("[\t ]+");
	
					if (n2id.containsKey(tokens[0]) && n2id.containsKey(tokens[1])) {
						i1 = n2id.get(tokens[0]);
						i2 = n2id.get(tokens[1]);
							
						double pw = Float.parseFloat(tokens[2]);
						if ((pw <= 0.05) && ((AMeGL[i1].size() > 0) && (AMeGL[i2].size() > 0))) {
							if (pw >= 1e-6)
								w = (- Math.log(pw) / Math.log(100.0));
							else
								w = 3.0;
	
							if (CoG[i1][i2] < 0.000001) {
								UG[i1].add(i2);
								UG[i2].add(i1);
								CoG[i1][i2] = (float) w;
								CoG[i2][i1] = CoG[i1][i2];
								ACoGL[i1].add(new Edge(i2, CoG[i1][i2]));
								ACoGL[i2].add(new Edge(i1, CoG[i1][i2]));
							}
						}
					}
				}
				br.close();
			}

			// read the graph input file again and fill out adjacent matrix and list	
			if (!funNetworkFile.isEmpty()){
				br = new BufferedReader(new InputStreamReader(new FileInputStream(funNetworkFile)));
	
				while ((line = br.readLine()) != null) {
					line = line.trim();
					String[] tokens = line.split("[\t ]");
	
					if (n2id.containsKey(tokens[0]) && n2id.containsKey(tokens[1])) {
						i1 = n2id.get(tokens[0]);
						i2 = n2id.get(tokens[1]);
						//w = Float.parseFloat(tokens[2]) / 1000.0f;
						w = Float.parseFloat(tokens[2]);
						
						if (G[i1][i2] < 0.000001) {
							UG[i1].add(i2);
							UG[i2].add(i1);
							G[i1][i2] = (float) w;
							G[i2][i1] = G[i1][i2];
							AGL[i1].add(new Edge(i2, G[i1][i2]));
							AGL[i2].add(new Edge(i1, G[i1][i2]));
						}
					}
				}
				br.close();
			}
		} catch (Exception e) {
			e.printStackTrace();

		}
		
		System.out.println("Succesfully read all the input files.\n");
		
		System.out.println("Starting to build Integer Linear Program model.\n");

		int count = 0;
		try {

			PrintWriter writer = new PrintWriter(outputLPFile, "UTF-8");
			writer.println("MAXIMIZE");
			writer.print(" obj: ");

			// float we;
			for (i = 0; i < n; i++)
				for (int j1 : UG[i]) {
					j = j1;

					if (method.equals("bemewithfun")) {
						if (i < j) {
							if ((MeG[i][j] > 0.0) || (G[i][j] > 0.0))
								writer.print(coeff(MeG[i][j] - G[i][j]) + " " + var("z", i, j));
						}
						if ((i < j) && (G[i][j] > 0.0))
							for (k = 1; k <= numClus - 1; k++)
								writer.print(coeff(G[i][j]) + " " + var("x", i, j, k));
					} else if (method.equals("bemewithco")) {
						if (i < j) {
							if ((MeG[i][j] > 0.0) || (CoG[i][j] > 0.0))
								writer.print(coeff(MeG[i][j] - CoG[i][j]) + " " + var("z", i, j));
						}
						if ((i < j) && ((MeG[i][j] > 0.0) || (CoG[i][j] > 0.0)))
							for (k = 1; k <= numClus - 1; k++)
								writer.print(coeff(CoG[i][j] - MeG[i][j]) + " " + var("x", i, j, k));
					} else if (method.equals("becowithmefun")) {
						if (i < j) {
							if ((MeG[i][j] > 0.0) || (CoG[i][j] > 0.0) || (G[i][j] > 0.0))
								writer.print(coeff(CoG[i][j] - MeG[i][j] - G[i][j]) + " " + var("z", i, j));
						}
						if ((i < j) && ((MeG[i][j] > 0.0) || (CoG[i][j] > 0.0) || (G[i][j] > 0.0)))
							for (k = 1; k <= numClus - 1; k++)
								writer.print(coeff(MeG[i][j] - CoG[i][j] + G[i][j]) + " " + var("x", i, j, k));
					} 

				}
			writer.println();
			writer.println("SUBJECT TO");


			count = 0;
			for (i = 0; i < n; i++) {
				for (int j1 : UG[i]) {
					j = j1;
					count++;
					writer.print(" ZijConstr" + count + ": ");
					writer.print(var("z", i, j));
					for (k = 1; k <= numClus - 1; k++)
						writer.print(" + " + var("x", i, j, k));
					writer.print(" + " + var("x", i, j, unselectedCluster));
					writer.print(" = 1 ");
					writer.println();
				}
			}

			count = 0;
			for (i = 0; i < n; i++) {
				for (int j1 : UG[i]) {
					j = j1;
					for (k = 1; k <= numClus - 1; k++) {
						count++;
						writer.print(" XijkConstr" + count + ": ");
						writer.print(var("x", i, j, k) + " - " + var("y", i, k) + " - " + var("y", j, k) + " >= -1");
						writer.println();
						count++;
						writer.print(" XijkConstr" + count + ": ");
						writer.print(var("x", i, j, k) + " - " + var("y", i, k) + " <= 0");
						writer.println();
						count++;
						writer.print(" XijkConstr" + count + ": ");
						writer.print(var("x", i, j, k) + " - " + var("y", j, k) + " <= 0");
						writer.println();
					}
					count++;
					writer.println(
							" XijkConstr" + count + ": " + var("x", i, j, unselectedCluster) + " - " + var("y", i, unselectedCluster) + " >= 0");
					count++;
					writer.println(
							" XijkConstr" + count + ": " + var("x", i, j, unselectedCluster) + " - " + var("y", j, unselectedCluster) + " >= 0");
					count++;
					writer.println(" XijkConstr" + count + ": " + var("x", i, j, unselectedCluster) + " - " + var("y", j, unselectedCluster) + " - "
							+ var("y", i, unselectedCluster) + " <= 0");
				}
			}

			//	the constraints ensure each gene is in one cluster
			count = 0;
			for (i = 0; i < n; i++) {
				count++;
				writer.print(" OneClusConstr" + count + ": ");
				writer.print(var("y", i, unselectedCluster));
				for (k = 1; k <= numClus - 1; k++)
					writer.print(" + " + var("y", i, k));
				writer.print(" = 1");
				writer.println();
			}

			//	each cluster has at most MAXMODULESIZE genes
			count = 0;
			for (k = 1; k <= numClus - 1; k++) {
				count++;
				writer.print(" SizeConstr" + count + ": ");
				writer.print(var("y", 0, k));
				for (i = 1; i < n; i++)
					writer.print(" + " + var("y", i, k));
				writer.print(" <= " + MAXMODULESIZE);
				writer.println();
			}

			//	each cluster has at least one gene
			for (k = 1; k <= numClus - 1; k++) {
				count++;
				writer.print(" SizeConstr" + count + ": ");
				writer.print(var("y", 0, k));
				for (i = 1; i < n; i++)
					writer.print(" + " + var("y", i, k));
				writer.println(" >= 1");
			}
			
			//	breaking symmetry constraint
			count=0;
			for (k = 1; k <= numClus - 1; k++) {
				count++;
				writer.print(" BreakSymConstr" + count + ": ");
				for (i = 0; i < n; i++)
					writer.print(" + " + i + " " +var("y", i, k));
				for (i = 0; i < n; i++)
					writer.print(" - " + i + " " +var("y", i, k+1));
				writer.println(" <= 0");
			}
			
			//	to enforce there are two other types of edges in the solution
			if (method.equals("becowithmefun")){
				writer.print(" becowithmefunConstr : ");
				for (i = 0; i < n; i++) 
					for (int j1 : UG[i]) {
						j = j1;
						if ((MeG[i][j] > 0.0)||(G[i][j] > 0.0))
						for (k = 1; k <= numClus - 1; k++) 
							writer.print(" + "+var("x", i, j, k));
					}
				writer.println(" >= 0.99");	
			}
			
			// to enforce the density of each module
			count = 0;
			float sum;
			for (k = 1; k <= numClus - 1; k++)
				for (i = 0; i < n; i++) {
					// if (AGL[i].size()>0){
					count++;
					double d = 0.699f;
					sum = (float) ((MAXMODULESIZE - 1.0f) * d);

					writer.print(" DensityConstr" + count + ": -" + (sum + d) + " " + var("y", i, k));

					HashSet<Integer> ai = new HashSet<Integer>();
					if ((method.equals("bemewithfun")) || (method.equals("becowithmefun")) || (method.equals("withmefun"))) {
						for (i1 = 0; i1 < AGL[i].size(); i1++)
							ai.add(AGL[i].get(i1).v);
					} else if (method.equals("bemewithco")) {
						for (i1 = 0; i1 < ACoGL[i].size(); i1++)
							ai.add(ACoGL[i].get(i1).v);
					} else if (method.equals("withme")) {
						for (i1 = 0; i1 < AMeGL[i].size(); i1++)
							ai.add(AMeGL[i].get(i1).v);
					}

					for (j = 0; j < n; j++)
						if (ai.contains(j))
							writer.print(" + " + (1.0 - d) + " " + var("y", j, k));
						else if (j != i)
							writer.print(" - " + d + " " + var("y", j, k));
					writer.println(" >= -" + (sum + d));
				}
			
			// to enforce for every gene in a module there is at least a between edge
			count = 0;
			for (i = 0; i < n; i++) {
				count++;
				writer.print(" ExAtLeastConstr" + count + ": ");
				for (k = 1; k <= numClus - 1; k++)
					writer.print(" - " + var("y", i, k));

				if (method.equals("becowithmefun")) {
					for (i1 = 0; i1 < ACoGL[i].size(); i1++) {
						j = ACoGL[i].get(i1).v;
						writer.print(" + " + var("z", i, j));
					}
				} else if ((method.equals("bemewithfun")) || (method.equals("bemewithco")))  {
					for (i1 = 0; i1 < AMeGL[i].size(); i1++) {
						j = AMeGL[i].get(i1).v;
						writer.print(" + " + var("z", i, j));
					}
				}
				writer.println(" >= 0");
			}
					
			writer.println("BOUNDS");

			for (i = 0; i < n; i++) 
				for (int j1 : UG[i]) {
					j = j1;
					for (k = 1; k <= numClus; k++)
						writer.println(" 0 <= " + var("x", i, j, k) + " <= 1");
				}
				
			for (i = 0; i < n; i++) 
				for (int j1 : UG[i]) {
					j = j1;
					writer.println(" 0 <= " + var("z", i, j) + " <= 1");
				}


			writer.println("BINARY");
			// indicator variable y_i_k	denoting the membership of gene i in cluster k   
			for (i = 0; i < n; i++) {
				for (k = 1; k <= numClus; k++)
					writer.println(var("y", i, k));
			}

			// indicator variable s_k_i denoting whether module k has i genes 
			for (k = 1; k <= numClus - 1; k++) {
				for (i = 1; i <= MAXMODULESIZE; i++) {
					writer.println(var("s", k, i));
				}
			}
			writer.println("END");
			writer.close();
			
			System.out.println("Finished building the ILP model.\n");
			
			String solutionLPFile=outputPrefix+"_"+(numClus-1)+"_modules.out";
			String cmd="";
			
			// removing all intermediate files
			ProcessBuilder pb;
			if (System.getProperty("os.name").startsWith("Win")){
				cmd="del "+outputPrefix+"_"+(numClus-1)+"_modules.out";
				pb = new ProcessBuilder("CMD", "/C", cmd);
			}
			else{
				cmd="rm "+outputPrefix+"_"+(numClus-1)+"_modules.out";
				pb = new ProcessBuilder("/bin/sh", "-c", cmd);
			}
            Process p = pb.start();
            p.waitFor();
            
            System.out.println("Calling CPLEX to solve ILP model.");
            System.out.println("Maximum solving time is two hours and all the threads will be used.");
            System.out.println("This might take a while...\n");
			
			// calling cplex to solve the output lp file
			String cplexCall="cplex";		
			cmd=cplexCall+" -c \"read " + outputLPFile +"\" \"set timelimit 7200\" \"set threads 0\" \"set mip ordertype 3\" \"mipopt\" \"write " + solutionLPFile + " sol\"";
			if (System.getProperty("os.name").startsWith("Win"))
				pb = new ProcessBuilder("CMD", "/C", cmd);
			else
				pb = new ProcessBuilder("/bin/sh", "-c", cmd);
            p = pb.start();
            p.waitFor();
            
            System.out.println("CPLEX finished solving ILP model.\n");
			
            System.out.println("Reading output soluton from CPLEX.\n");
    		count=0;
    		int clus;
    		String outputResultFile=outputPrefix+"_"+(numClus-1)+"_modules";
    		HashSet<Integer>[] modules=new HashSet[numClus-1];
    		HashSet<Integer> selected=new HashSet<Integer>();
    		HashMap<Integer,Integer> hm=new HashMap<Integer,Integer>();
    		for (i=0;i<numClus-1;i++)
    			modules[i]=new HashSet<Integer>();
    		
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(solutionLPFile)));
			String pattern="y_([0-9]*)_([0-9]*)";
			String yvalues="value=\"(.*?)\"/>";
			String ovalues="objectiveValue=\"(.*?)\"";
			double ofv=0.0;
			while ((line = br.readLine()) != null) {
				Matcher m = Pattern.compile(pattern).matcher(line);
				Matcher m1 = Pattern.compile(yvalues).matcher(line);
				Matcher m2 = Pattern.compile(ovalues).matcher(line);
				
				if (m2.find()){
					ofv=Double.parseDouble(m2.group(1));
				}
				if (m.find()&&m1.find()) {
					String[] tokens = m.group(0).split("_");
					if (Double.parseDouble(m1.group(1))>0.99){
					clus=Integer.parseInt(tokens[2]);
						if ((clus<numClus)&&(clus>0)){
							modules[clus-1].add(Integer.parseInt(tokens[1]));
							selected.add(Integer.parseInt(tokens[1]));
							hm.put(Integer.parseInt(tokens[1]), clus-1);
						}
					}
				}
			}
			br.close();
			
			System.out.println("Writing summary of the result to "+outputResultFile+".\n");
			
			writer = new PrintWriter(outputResultFile, "UTF-8");
			
			// print out the objective function value
			writer.println("Score = "+ofv);
			writer.println();
			
			// print out the genes in each module
			for (k=0;k<numClus-1;k++){
				writer.print("Module "+(k+1)+":");
				for (Integer u1:modules[k]){
					writer.print(" "+id2n.get(u1));
				}
				writer.println();
			}
			writer.println();
			
			//print out within-module edges
			writer.println("Within-module edges:");
			for (Integer u1:selected)
			for (Integer v1:selected)
			if ((u1<v1)&&(hm.get(u1)==hm.get(v1))){
				if ((G[u1][v1]>0.0)&&((method.equals("bemewithfun"))||(method.equals("becowithmefun"))))	
					writer.println(id2n.get(u1)+","+id2n.get(v1)+",pp,"+1.0);
				if ((MeG[u1][v1]>0.0)&&(method.equals("becowithmefun")))
					writer.println(id2n.get(u1)+","+id2n.get(v1)+",me,"+MeG[u1][v1]);
				if ((CoG[u1][v1]>0.0)&&(method.equals("bemewithco")))
					writer.println(id2n.get(u1)+","+id2n.get(v1)+",co,"+CoG[u1][v1]);
			}
			writer.println();
			
			//	print out between-module edges
			writer.println("Between-module edges:");
			for (Integer u1:selected)
			for (Integer v1:selected)
			if ((u1<v1)&&(hm.get(u1)!=hm.get(v1))){
				if ((MeG[u1][v1]>0.0)&&((method.equals("bemewithfun"))||(method.equals("bemewithco"))))	
					writer.println(id2n.get(u1)+","+id2n.get(v1)+",me,"+MeG[u1][v1]);
				if ((CoG[u1][v1]>0.0)&&(method.equals("becowithmefun")))
					writer.println(id2n.get(u1)+","+id2n.get(v1)+",co,"+CoG[u1][v1]);
			}
			
			writer.close();		
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

