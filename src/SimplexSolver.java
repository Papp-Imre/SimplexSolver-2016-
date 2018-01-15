import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

public class SimplexSolver {
	private int N = 0;
	private int M = 0;
	
	private float A[][];
	
	private float originalArray[][];
	private float tempArray[][];
	
	private float c[];
	private float tempC[];
	
	private float b[];
	private float tempB[];
	
	private int u[];
	private int x[];
	
	private float genValue = 0;
	private float newGenValue = 0;
	
	private float result = 0;
	private String filePath = "szimplexinput.txt";
	
	private float[] originalB;
	private float[] originalC;
    
    public static void main(String[] args) {
    	
    	SimplexSolver solver = new SimplexSolver();
    	if(args.length > 0) {
    		solver.filePath = args[0];
    	}
        solver.start();
    }
    
    private void start() {
    	loadFile();
        solve();
    }
    
    private void solve() {
    	createTempArrays();
    	
    	int colIndex = -2; // init
    	
    	int num = 0;
    	while(true) {
    		try {
				colIndex = getMaxIndexC();
			} catch (Exception e) {
				System.out.println("M�r nincs olyan c elem, ami felett v�laszthat� lenne gener�l� elem, viszont maradt m�g pozit�v c, �gy nincs megold�s.");
				break;
			}
    		
    		System.out.println("------------------------");

    		if(colIndex == -1) {
    			// ide jutunk ha minden C elem negat�v
    			// hibas c input - TODO: consolera
    			break;
    		}
    		
    		System.out.println("");
    		System.out.println("Gener�l� elem oszlopa: "+(colIndex+1));
    		
    		// szukkeresztmetszet
    		int rowIndex = getRowResult(colIndex);

    		System.out.println("Gener�l� elem sora: "+(rowIndex+1));
    		System.out.println("");
    		
    		// generalo elem szamolasa
    		calcGenValue(rowIndex, colIndex);
    		
    		// generalo elem soranak szamolasa
    		calcNewRow(rowIndex, colIndex);
    		
    		// generalo elem oszlopanak szamolasa
    		calcNewColumn(rowIndex, colIndex);
    		
    		// matrix tobbi elemenek szamolasa
    		calcNewMatrix(rowIndex, colIndex);

    		// b vektor sz�mol�sa
    		calcNewB(rowIndex, colIndex);

    		// c vektor sz�mol�sa
    		calcNewC(rowIndex, colIndex);

    		// c�lf�ggv�ny �rt�k�nek sz�mol�sa
    		calcNewResult(rowIndex, colIndex);

    		// az �rhat� t�mb�k elment�se az olvashat�akba
    		saveTempArray();
    		
    		num++;
    		System.out.println("A "+num+". l�p�s. Gener�l� elem oszlopa: "+(colIndex+1));
    		System.out.println("A "+num+". l�p�s. Gener�l� elem sora: "+(rowIndex+1));
    		System.out.println("A "+num+". l�p�s. Gener�l� elem �rt�ke: "+A[rowIndex][colIndex]);
    		System.out.println("");
    		
    		// Print matrix
    		System.out.println("Print matrix: ");
    		for (int i=0; i < A.length; i++){
    			printArray(A[i]);
    		}
    		
    		// Print b
    		System.out.println("");
    		System.out.println("Print b: ");
    		printArray(b);
    		
    		// Print c
    		System.out.println("");
    		System.out.println("Print c: ");
    		printArray(c);
    		
    		// Print u
    		System.out.println("");
    		System.out.println("Print u: ");
    		printArray(u);
    		
    		// Print x
    		System.out.println("");
    		System.out.println("Print x: ");
    		printArray(x);
    		
    		// Print result
    		System.out.println("");
    		System.out.println("Print result: "+result);
    		System.out.println("");
    	}
    	System.out.println("");
    	System.out.println("Megold�s megtal�lva!");
    	System.out.println("");
    	System.out.println("V�geredm�ny: "+(result*-1));
    	System.out.println("");

    	float[] xArray = new float[x.length];
    	float[] uArray = new float[u.length];
    	
    	for (int i = 0; i < u.length; i++) {
    		int idx = u[i];
    		float value = 0;
			if(u[i]<0) {
				idx = (idx * -1) - 1;
				value = b[i];
				xArray[idx] = value;
			} else {
				uArray[i] = value;
			}
		}
    	
    	for (int i = 0; i < x.length; i++) {
    		int idx = x[i];
    		float value = 0;
			if(x[i]<0) {
				idx = (idx * -1) - 1;
				value = (c[i]*-1);
				uArray[idx] = value;
			} 
		}
    	
    	for (int i = 0; i < xArray.length; i++) {
			System.out.println("Az x"+(i+1)+" elem �rt�ke: "+xArray[i]);
		}
    	System.out.println("");    	
    	for (int i = 0; i < uArray.length; i++) {
			System.out.println("Az u"+(i+1)+" elem �rt�ke: "+uArray[i]);
    	}
    	System.out.println("");
    	for (int i = 0; i < u.length; i++) {
    		float sum = 0;
    		for (int j = 0; j < x.length; j++) {
				sum+= xArray[j]*originalArray[i][j];
			}
    		System.out.println("A b"+(i+1)+"-b�l kihaszn�lt kapacit�s: "+sum+", megmaradt: "+ (originalB[i]-sum));
		}
    	System.out.println("");
    	calculateAlternativeResult();
    }

	private void calculateAlternativeResult() {
		System.out.println("----------------------------");
		System.out.println("");
		System.out.println("Alternat�v optimum keres�se:");
		System.out.println("");
		int colIndex = -1;
		try {
			colIndex = getMaxIndexC(true);
		} catch (Exception e) {
			System.out.println("");
			System.out.println("Nincs alternat�v megold�s.");
			return;
		}
		
		if(colIndex == -1) {
			System.out.println("");
			System.out.println("Nincs alternat�v megold�s.");
			return;
		}

		System.out.println("");
		System.out.println("Gener�l� elem oszlopa: "+(colIndex+1));
		
		// sz�k keresztmetszet sz�mol�sa
		int rowIndex = getRowResult(colIndex);

		System.out.println("Gener�l� elem sora: "+(rowIndex+1));
		System.out.println("");
		
		// gener�l� elem sz�mol�sa
		calcGenValue(rowIndex, colIndex);
		
		// gener�l� elem sor�nak sz�mol�sa
		calcNewRow(rowIndex, colIndex);
		
		// gener�l� elem oszlop�nak sz�mol�sa
		calcNewColumn(rowIndex, colIndex);
		
		// m�trix t�bbi elem�nek sz�mol�sa
		calcNewMatrix(rowIndex, colIndex);

		// b vektor sz�mol�sa
		calcNewB(rowIndex, colIndex);

		// c vektor sz�mol�sa
		calcNewC(rowIndex, colIndex);

		// c�lf�ggv�ny �rt�k�nek sz�mol�sa
		calcNewResult(rowIndex, colIndex);

		// az �rhat� t�mb�k elment�se az olvashat�akba
		saveTempArray();
		
		System.out.println("A alternat�v optimum keres�se k�zben. Gener�l� elem oszlopa: "+(colIndex+1));
		System.out.println("A alternat�v optimum keres�se k�zben. Gener�l� elem sora: "+(rowIndex+1));
		System.out.println("A alternat�v optimum keres�se k�zben. Gener�l� elem �rt�ke: "+A[rowIndex][colIndex]);
		System.out.println("");
		System.out.println("Alternat�v optimum megtal�lva!");
		System.out.println("");
		System.out.println("V�geredm�ny: "+(result*-1));
    	System.out.println("");

		// x, u vektor �rt�keinek lement�se t�mb�kbe
    	float[] xArray = new float[x.length];
    	float[] uArray = new float[u.length];
    	
    	for (int i = 0; i < u.length; i++) {
    		int idx = u[i];
    		float value = 0;
			if(u[i]<0) {
				idx = (idx * -1) - 1;
				value = b[i];
				xArray[idx] = value;
			} else {
				uArray[i] = value;
			}
		}
    	
    	for (int i = 0; i < x.length; i++) {
    		int idx = x[i];
    		float value = 0;
			if(x[i]<0) {
				idx = (idx * -1) - 1;
				value = (c[i]*-1);
				uArray[idx] = value;
			}
		}
    	
    	// ki�rat�sok, x,u vektorok �rt�kei, marad�ksz�ml�l�s
    	for (int i = 0; i < xArray.length; i++) {
			System.out.println("Az x"+(i+1)+" elem �rt�ke: "+xArray[i]);
		}
    	System.out.println("");
    	for (int i = 0; i < uArray.length; i++) {
			System.out.println("Az u"+(i+1)+" elem �rt�ke: "+uArray[i]);
    	}
    	System.out.println("");
    	for (int i = 0; i < u.length; i++) {
    		float sum = 0;
    		for (int j = 0; j < x.length; j++) {
				sum+= xArray[j]*originalArray[i][j];
			}
    		System.out.println("A b"+(i+1)+"-b�l kihaszn�lt kapacit�s: "+sum+", megmaradt: "+ (originalB[i]-sum));
		}
    	System.out.println("");
	}

	/**
     * Kisz�molja a maxim�lisan v�laszthat� c �rt�ket 
     * @return c elem indexe
     * @throws Exception - ha nincs v�laszthat� elem, akkor hib�t dob
     */
	private int getMaxIndexC() throws Exception {
		return getMaxIndexC(false);
	}

	/**
     * Kisz�molja a maxim�lisan v�laszthat� c �rt�ket 
     * @param zeroEnabled - a 0 megengedett c elem-e a v�laszt�shoz 
     * @return c elem indexe
     * @throws Exception - ha nincs v�laszthat� elem, akkor hib�t dob
	 */
	private int getMaxIndexC(boolean zeroEnabled) throws Exception{

		int colIndex = -1;
		
    	int[] possibleCIndexes = getIndexesBySortC(zeroEnabled);
    	
    	int actualIndexNumber = 0;
    	if(possibleCIndexes.length > 0) {
    		boolean isValidC = false; 
    		
    		while(!isValidC) {
    			colIndex = possibleCIndexes[actualIndexNumber];
    			actualIndexNumber++;
    			
    			isValidC = checkCIsValid(colIndex);
    			
    			if(!isValidC) {
    				System.out.println("Az �sszes "+(colIndex+1)+". c elem feletti sorelem negat�v, de nem lehet negat�v elemet gener�l� elemnek v�lasztani, �gy az a c elem nem v�laszthat�.");
    			}
    		}
    	}
    	
    	return colIndex;
	}

	/**
	 * Megvizsg�lja, hogy az adott index� (eset�nkben a gener�l� elem oszlop�ban l�v�) sz�mok �rv�nyesek, teh�t nemnegat�vak-e. 
	 * @param colIndex - gener�l� elem oszlop�nak sz�ma
	 * @return �rv�nyes-e az oszlop sz�k keresztmetszethez
	 */
	private boolean checkCIsValid(int colIndex) {
		float[] column = getColumn(colIndex, tempArray);
		int notValidNumbers = 0;
		
		for(int i = 0; i < column.length; i++) {
			if(column[i] <= 0) {
				notValidNumbers++;
			}
		}
		
		return (notValidNumbers != column.length);
	}

	/**
	 * V�ltoz�sok ny�lv�ntart�s�hoz az �rhat� t�mb�k l�trehoz�sa
	 */
	private void createTempArrays() {
		tempArray = deepCopy(A);
		tempB = Arrays.copyOf(b, b.length);
		tempC = Arrays.copyOf(c, c.length);
	}
	
	/**
	 * V�ltoz�sok elment�se
	 */
	private void saveTempArray() {
		A = deepCopy(tempArray);
		b = Arrays.copyOf(tempB, tempB.length);
		c = Arrays.copyOf(tempC, tempC.length);
	}

	/**
	 * Kisz�molja a gener�l� elemet, illetve elv�gzi az x �s u vektorok cser�j�t
	 * @param rowIndex
	 * @param colIndex
	 */
	private void calcGenValue(int rowIndex, int colIndex) {
		genValue = tempArray[rowIndex][colIndex];
		tempArray[rowIndex][colIndex]= 1/(tempArray[rowIndex][colIndex]);
		newGenValue = tempArray[rowIndex][colIndex];
		
		int tmp = x[colIndex];
		x[colIndex] = (u[rowIndex]*-1);
		u[rowIndex] = (tmp*-1); 
		
		System.out.println("A gener�l� elem eredetileg: "+genValue+" �s m�dosulva: "+newGenValue);
		System.out.println("");
	}
	
	/**
	 * Kisz�molja a gener�l� elem oszlop�ban l�v� sz�mokat a m�trixban.
	 * @param genIndex
	 * @param colIndex
	 */
	private void calcNewColumn(int genIndex, int colIndex) {
		float[] column = getColumn(colIndex, tempArray);
		float[] aColumn = getColumn(colIndex, A);
		 for(int i = 0; i < column.length; i++) {
			 if (column[i] == aColumn[i] && i!=genIndex) {
				 column[i] = column[i]*(-1*column[genIndex]);
//				 tempArray[i][colIndex] = column[i];
				 if(column[i]==-0.0f) {
//				 tempArray[i][colIndex] = 0;
					 column[i] = 0;
				 }
				 System.out.println("M�dosult oszlopelemek �rt�ke: "+(i+1)+". elemb�l lett: "+column[i]);
			 }

		 }
		 for(int k=0; k < column.length; k++) {
			 tempArray[k][colIndex] = column[k];
		 }
	}

	/**
	 * Kisz�molja a gener�l� elem sor�ban l�v� sz�mokat a m�trixban.
	 * @param rowIndex
	 * @param genIndex
	 */
	private void calcNewRow(int rowIndex, int genIndex) {
		float[] row = getRow(tempArray, rowIndex);
		float[] aRow = getRow(A, rowIndex);
		 for(int i = 0; i < row.length; i++) {
			 if (row[i] == aRow[i]) {
				 row[i] = row[i]*row[genIndex];
				 System.out.println("M�dosult sorelemek �rt�ke: "+(i+1)+". elemb�l lett: "+row[i]);
			 }
			 if(row[i]==-0.0f) {
				 row[i] = 0;
			 }
			 
		 }
		 tempArray[rowIndex] = row;
		 System.out.println("");
	}
	
	/**
	 * Kisz�molja a b vektort.
	 * @param rowIndex
	 * @param columnIndex
	 */
	private void calcNewB(int rowIndex, int columnIndex) {
		float[] row = tempB;
		 for(int i = 0; i < row.length; i++) {
			 if (i == rowIndex) {
				 row[i] = row[i] * newGenValue;
				 System.out.println("M�dosult B (gen.elem sor�ban) elem �rt�ke: "+(i+1)+". elemb�l lett: "+row[i]);
			 }
			 else{
				 row[i] = row[i]-(-b[rowIndex]*tempArray[i][columnIndex]);
				 System.out.println("M�dosult B elemek �rt�ke: "+(i+1)+". elemb�l lett: "+row[i]);
			 }
			 if(row[i]==-0.0f) {
				 row[i] = 0;
			 }
		 }
		 tempB = row;
		 System.out.println("");
	}

	/**
	 * Kisz�molja a c vektort.
	 * @param rowIndex
	 * @param columnIndex
	 */
	private void calcNewC(int rowIndex, int columnIndex) {
		float[] col = tempC;
		 for(int i = 0; i < col.length; i++) {
			 if (i == columnIndex) {
				 col[i] = col[i] * (-1 * newGenValue);
				 System.out.println("M�dosult C (gen.elem oszlop�ban) elem �rt�ke: "+(i+1)+". elemb�l lett: "+col[i]);
			 }
			 else{
				 col[i] = col[i]-(c[columnIndex]*tempArray[rowIndex][i]);
				 System.out.println("M�dosult C elemek �rt�ke: "+(i+1)+". elemb�l lett: "+col[i]);
			 }
			 if(col[i]==-0.0f) {
					col[i] = 0;
			 }
		 }
		 tempC = col;
		 System.out.println("");
	}
	
	/**
	 * Kisz�molja a c�lf�ggv�ny �rt�k�t
	 * @param rowIndex
	 * @param colIndex
	 */
	private void calcNewResult(int rowIndex, int colIndex) {
		float[] col = tempB;
		float[] row = c;
		result = result - row[colIndex] * col[rowIndex];
		System.out.println("C�lf�ggv�ny �j �rt�ke: "+result);
		System.out.println("");
	}
		 
	/**
	 * Kisz�molja a m�trix t�bbi elem�t.
	 * @param rowIndex
	 * @param colIndex
	 */
	private void calcNewMatrix(int rowIndex, int colIndex) {
		float[][] oldMatrix = A;
		float[][] newMatrix = tempArray;
		for (int i=0; i < newMatrix.length; i++){
			for (int j=0; j < newMatrix[i].length; j++){
				// a kor�bban kit�lt�tt (gener�l� elem sora, oszlopa) mez�k figyelmen k�v�l hagy�sa
				if (i == rowIndex || j == colIndex){
					continue;
				}
				// �j koordin�ta = r�gi koordin�ta - (gener�l� elem sor�nak �j �rt�ke * gener�l� elem r�gi oszlop�b�l a megfelel� koordin�ta)
				newMatrix[i][j] = oldMatrix[i][j] - (newMatrix[rowIndex][j] * oldMatrix[i][colIndex]);
				if(newMatrix[i][j]==-0.0f) {
					newMatrix[i][j] = 0;
				}
			}
		}
		// ment�s
		tempArray = newMatrix;
		System.out.println("");
	}

	/**
	 * F�jlb�l beolvas�s, a mintak�nt kapott input alapj�n dinamikusan (n * m)
	 */
	private void loadFile() {
		FileInputStream fstream = null;
		BufferedReader br = null;
		try{
			// filePath - alapvet�en a fent megadott el�r�si �t, futtat�si param�ter megad�sa eset�n pedig az
            fstream = new FileInputStream(filePath );
            br = new BufferedReader(new InputStreamReader(fstream));
            String strLine;
            N = Integer.parseInt(br.readLine());
            M = Integer.parseInt(br.readLine());
            br.readLine();
            
            // megfelel� m�ret� t�mb�k l�trehoz�sa a t�rol�shoz
            initArrays();
            
            int i = 0;
            int j = 0;
            for(int k = 0; k < M; k++) {
                strLine = br.readLine();
                StringTokenizer st = new StringTokenizer(strLine," ");
                j = 0;
                while(st.hasMoreTokens()){                    
                	String token = st.nextToken();
                	A[i][j] = Float.parseFloat(token);
                	j++;
                }
                i++;
            }
            
            br.readLine();
            
            strLine = br.readLine();
            
            StringTokenizer st =new StringTokenizer(strLine, " ");
            int index = 0;
            while (st.hasMoreTokens()) {
				String token = st.nextToken();
				c[index] = Float.parseFloat(token);
				index++;
            }
            
            br.readLine();
            
            strLine = br.readLine();
            index = 0;
            st =new StringTokenizer(strLine, " ");
            while (st.hasMoreTokens()) {
				String token = st.nextToken();
				b[index] = Float.parseFloat(token);
				index++;
            }
            
            br.close();
            
            // eredeti m�trix, vektorok elment�se (deep copy)
            originalArray = deepCopy(A);
            originalB = Arrays.copyOf(b, b.length);
            originalC = Arrays.copyOf(c, c.length);
    	    
		
		}		
        catch(IOException ioe){
            System.err.println(ioe.getMessage());
        } finally {
			if(fstream != null) {
				try {
					fstream.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			if(br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Megfelel� m�ret� t�mb�k l�trehoz�sa a t�rol�shoz.
	 */
	private void initArrays() {
		A = new float[M][N];
		
		c = new float[N];
		
		b = new float[M];
		
		u = new int[M];
		
		x = new int[N];
		
		for (int i = 1; i <= u.length; i++) {
			u[i-1] = i;
		}
		
		for (int i = 1; i <= x.length; i++) {
			x[i-1] = i;
		}
	}
	
	/**
	 * K�tdimenzi�s t�mb teljes m�sol�sa egy �jba.
	 * @param original - az eredeti 2d t�mb
	 * @return a m�solat
	 */
	public static float[][] deepCopy(float[][] original) {
	    if (original == null) {
	        return null;
	    }

	    final float[][] result = new float[original.length][];
	    for (int i = 0; i < original.length; i++) {
	        result[i] = Arrays.copyOf(original[i], original[i].length);
	    }
	    return result;
	}
    
	/**
	 * Az A m�trix ki�r�sa a konzolra.
	 */
	private void printMatrix() {
		for(int i = 0; i < A.length; i++) {
			for(int j = 0; j < A[i].length; j++) {
				System.out.print(A[i][j]+" ");
			}
			System.out.println("");
		}
	}
	
	/**
	 * Param�terben kapott vektor ki�r�sa a konzolra.
	 * @param a - eg�sz sz�mokb�l �ll� vektor (u, v)
	 */
	private void printArray(int[] a) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i]+" ");
		}
		System.out.println("");
	}

	/**
	 * Param�terben kapott vektor ki�r�sa a konzolra.
	 * @param a - lebeg�pontos sz�mokb�l �ll� vektor (pl. b, c)
	 */
	private void printArray(float a[]) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i]+" ");
		}
		System.out.println("");
	}

	/**
	 * Visszaadja a pozit�v c elemek list�j�t.
	 * Rendezi a c vektort, majd a v�laszthat� �rt�keket bele�rj�k egy t�mbbe.
	 * @return - index t�mb
	 */
	public int[] getIndexesBySortC(boolean zeroEnabled) {
		float[] tempC = Arrays.copyOf(c, c.length);
		
		float helper = 0;
		// tempC sorbarendezese
		for (int i = 0; i < tempC.length; i++)
        {
            for (int j = 0; j < tempC.length; j++)
            {
                if (tempC[i] > tempC[j])
                {
                    helper = tempC[i];
                    tempC[i] = tempC[j];
                    tempC[j] = helper;
                }
            }
        }
		
		// negativ elemek keresese
		List<Integer> needDelete = new ArrayList<Integer>();
		for(int i = 0; i < tempC.length; i++) {
			if(zeroEnabled) {
				if(tempC[i] < 0) {
					needDelete.add(i);
					System.out.println("A "+tempC[i]+" < 0, ez�rt biztosan nem v�laszthat� a c - b�l");
				}
			} else {
				if(tempC[i] <= 0) {
					needDelete.add(i);
					System.out.println("A "+tempC[i]+" <= 0, ez�rt biztosan nem v�laszthat� a c - b�l");
				}
			}
		}
		
		int tmpSize = (tempC.length-needDelete.size());
		float[] tmp = new float[tmpSize];
		int tmpIdx = 0;
		
		//negativ elemek "torlese"
		for(int i = 0; i < tempC.length; i++) {
			if(needDelete.contains(i)) {
				i++;
			} else {
				tmp[tmpIdx] = tempC[i];
				tmpIdx++;
			}
		}
		
		tempC = tmp;
		int[] indexArray = new int[tempC.length];
		
		// indexek megkeresese	az eredeti c-ben
		if(tempC.length < 1) return indexArray;
		float value = tempC[0];
		int helperIndex = 1;
		for(int i = 0; i < tempC.length; i++) {
			if(tempC[i] != value) {
				value = tempC[i];
				helperIndex = 0;
			}
			
			for(int j = helperIndex; j < c.length; j++) {
				if(c[j] == tempC[i]) {
					indexArray[i] = j;
					helperIndex = j + 1;
					break;
				}
			}
		}
		
		return indexArray;
	}
	
	/**
	 * Sz�k keresztmetszet sz�m�t�sa 
	 * @param matrixIndex
	 * @return - legpozit�vabb v�laszt�s indexe
	 */
	public int getRowResult(int matrixIndex) {
		float[] column = getColumn(matrixIndex, tempArray);
		float max = getDivisonResultForRow(column[0], 0);
		int maxIndex = 0;
		float result = -1;
		for(int i = 1; i < column.length; i++) {
			if(column[i]<1) {
				continue;
			}
			result = getDivisonResultForRow(column[i], i);
			if(max > result || max < 0) {
				max = result;
				maxIndex = i;
			} 
		}
		
		return maxIndex;
	}
	
	

	/**
	 * K�sz� egy m�solatot a kapott 2dimenzi�s t�mb k�v�nt oszlop�r�l.
	 * @param matrixIndex - a k�rt oszlop indexe
	 * @param source - 2 dimenzi�s t�mb
	 * @return
	 */
	private float[] getColumn(int matrixIndex, float[][] source) {
		float[] retArray = new float[M];
		for(int i = 0; i < M; i++) {
			retArray[i] = source[i][matrixIndex];
	}
		return retArray;
	}
	
	/**
	 * Visszaadja a 2 dimenzi�s t�mb k�v�nt sor�t.
	 * @param source - 2 dimenzi�s t�mb
	 * @param matrixIndex - a k�rt sor indexe
	 * @return
	 */
	private float[] getRow(float[][] source, int matrixIndex) {
		return source[matrixIndex];
	}

	/**
	 * Sz�k keresztmetszethez sz�ks�ges oszt�s
	 * @param value
	 * @param index
	 * @return
	 */
	public float getDivisonResultForRow(float value, int index) {
		return b[index]/value;
	}
}
