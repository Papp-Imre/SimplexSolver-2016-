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
				System.out.println("Már nincs olyan c elem, ami felett választható lenne generáló elem, viszont maradt még pozitív c, így nincs megoldás.");
				break;
			}
    		
    		System.out.println("------------------------");

    		if(colIndex == -1) {
    			// ide jutunk ha minden C elem negatív
    			// hibas c input - TODO: consolera
    			break;
    		}
    		
    		System.out.println("");
    		System.out.println("Generáló elem oszlopa: "+(colIndex+1));
    		
    		// szukkeresztmetszet
    		int rowIndex = getRowResult(colIndex);

    		System.out.println("Generáló elem sora: "+(rowIndex+1));
    		System.out.println("");
    		
    		// generalo elem szamolasa
    		calcGenValue(rowIndex, colIndex);
    		
    		// generalo elem soranak szamolasa
    		calcNewRow(rowIndex, colIndex);
    		
    		// generalo elem oszlopanak szamolasa
    		calcNewColumn(rowIndex, colIndex);
    		
    		// matrix tobbi elemenek szamolasa
    		calcNewMatrix(rowIndex, colIndex);

    		// b vektor számolása
    		calcNewB(rowIndex, colIndex);

    		// c vektor számolása
    		calcNewC(rowIndex, colIndex);

    		// célfüggvény értékének számolása
    		calcNewResult(rowIndex, colIndex);

    		// az írható tömbök elmentése az olvashatóakba
    		saveTempArray();
    		
    		num++;
    		System.out.println("A "+num+". lépés. Generáló elem oszlopa: "+(colIndex+1));
    		System.out.println("A "+num+". lépés. Generáló elem sora: "+(rowIndex+1));
    		System.out.println("A "+num+". lépés. Generáló elem értéke: "+A[rowIndex][colIndex]);
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
    	System.out.println("Megoldás megtalálva!");
    	System.out.println("");
    	System.out.println("Végeredmény: "+(result*-1));
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
			System.out.println("Az x"+(i+1)+" elem értéke: "+xArray[i]);
		}
    	System.out.println("");    	
    	for (int i = 0; i < uArray.length; i++) {
			System.out.println("Az u"+(i+1)+" elem értéke: "+uArray[i]);
    	}
    	System.out.println("");
    	for (int i = 0; i < u.length; i++) {
    		float sum = 0;
    		for (int j = 0; j < x.length; j++) {
				sum+= xArray[j]*originalArray[i][j];
			}
    		System.out.println("A b"+(i+1)+"-bõl kihasznált kapacitás: "+sum+", megmaradt: "+ (originalB[i]-sum));
		}
    	System.out.println("");
    	calculateAlternativeResult();
    }

	private void calculateAlternativeResult() {
		System.out.println("----------------------------");
		System.out.println("");
		System.out.println("Alternatív optimum keresése:");
		System.out.println("");
		int colIndex = -1;
		try {
			colIndex = getMaxIndexC(true);
		} catch (Exception e) {
			System.out.println("");
			System.out.println("Nincs alternatív megoldás.");
			return;
		}
		
		if(colIndex == -1) {
			System.out.println("");
			System.out.println("Nincs alternatív megoldás.");
			return;
		}

		System.out.println("");
		System.out.println("Generáló elem oszlopa: "+(colIndex+1));
		
		// szûk keresztmetszet számolása
		int rowIndex = getRowResult(colIndex);

		System.out.println("Generáló elem sora: "+(rowIndex+1));
		System.out.println("");
		
		// generáló elem számolása
		calcGenValue(rowIndex, colIndex);
		
		// generáló elem sorának számolása
		calcNewRow(rowIndex, colIndex);
		
		// generáló elem oszlopának számolása
		calcNewColumn(rowIndex, colIndex);
		
		// mátrix többi elemének számolása
		calcNewMatrix(rowIndex, colIndex);

		// b vektor számolása
		calcNewB(rowIndex, colIndex);

		// c vektor számolása
		calcNewC(rowIndex, colIndex);

		// célfüggvény értékének számolása
		calcNewResult(rowIndex, colIndex);

		// az írható tömbök elmentése az olvashatóakba
		saveTempArray();
		
		System.out.println("A alternatív optimum keresése közben. Generáló elem oszlopa: "+(colIndex+1));
		System.out.println("A alternatív optimum keresése közben. Generáló elem sora: "+(rowIndex+1));
		System.out.println("A alternatív optimum keresése közben. Generáló elem értéke: "+A[rowIndex][colIndex]);
		System.out.println("");
		System.out.println("Alternatív optimum megtalálva!");
		System.out.println("");
		System.out.println("Végeredmény: "+(result*-1));
    	System.out.println("");

		// x, u vektor értékeinek lementése tömbökbe
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
    	
    	// kiíratások, x,u vektorok értékei, maradékszámlálás
    	for (int i = 0; i < xArray.length; i++) {
			System.out.println("Az x"+(i+1)+" elem értéke: "+xArray[i]);
		}
    	System.out.println("");
    	for (int i = 0; i < uArray.length; i++) {
			System.out.println("Az u"+(i+1)+" elem értéke: "+uArray[i]);
    	}
    	System.out.println("");
    	for (int i = 0; i < u.length; i++) {
    		float sum = 0;
    		for (int j = 0; j < x.length; j++) {
				sum+= xArray[j]*originalArray[i][j];
			}
    		System.out.println("A b"+(i+1)+"-bõl kihasznált kapacitás: "+sum+", megmaradt: "+ (originalB[i]-sum));
		}
    	System.out.println("");
	}

	/**
     * Kiszámolja a maximálisan választható c értéket 
     * @return c elem indexe
     * @throws Exception - ha nincs választható elem, akkor hibát dob
     */
	private int getMaxIndexC() throws Exception {
		return getMaxIndexC(false);
	}

	/**
     * Kiszámolja a maximálisan választható c értéket 
     * @param zeroEnabled - a 0 megengedett c elem-e a választáshoz 
     * @return c elem indexe
     * @throws Exception - ha nincs választható elem, akkor hibát dob
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
    				System.out.println("Az összes "+(colIndex+1)+". c elem feletti sorelem negatív, de nem lehet negatív elemet generáló elemnek választani, így az a c elem nem választható.");
    			}
    		}
    	}
    	
    	return colIndex;
	}

	/**
	 * Megvizsgálja, hogy az adott indexû (esetünkben a generáló elem oszlopában lévõ) számok érvényesek, tehát nemnegatívak-e. 
	 * @param colIndex - generáló elem oszlopának száma
	 * @return érvényes-e az oszlop szûk keresztmetszethez
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
	 * Változások nyílvántartásához az írható tömbök létrehozása
	 */
	private void createTempArrays() {
		tempArray = deepCopy(A);
		tempB = Arrays.copyOf(b, b.length);
		tempC = Arrays.copyOf(c, c.length);
	}
	
	/**
	 * Változások elmentése
	 */
	private void saveTempArray() {
		A = deepCopy(tempArray);
		b = Arrays.copyOf(tempB, tempB.length);
		c = Arrays.copyOf(tempC, tempC.length);
	}

	/**
	 * Kiszámolja a generáló elemet, illetve elvégzi az x és u vektorok cseréjét
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
		
		System.out.println("A generáló elem eredetileg: "+genValue+" és módosulva: "+newGenValue);
		System.out.println("");
	}
	
	/**
	 * Kiszámolja a generáló elem oszlopában lévõ számokat a mátrixban.
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
				 System.out.println("Módosult oszlopelemek értéke: "+(i+1)+". elembõl lett: "+column[i]);
			 }

		 }
		 for(int k=0; k < column.length; k++) {
			 tempArray[k][colIndex] = column[k];
		 }
	}

	/**
	 * Kiszámolja a generáló elem sorában lévõ számokat a mátrixban.
	 * @param rowIndex
	 * @param genIndex
	 */
	private void calcNewRow(int rowIndex, int genIndex) {
		float[] row = getRow(tempArray, rowIndex);
		float[] aRow = getRow(A, rowIndex);
		 for(int i = 0; i < row.length; i++) {
			 if (row[i] == aRow[i]) {
				 row[i] = row[i]*row[genIndex];
				 System.out.println("Módosult sorelemek értéke: "+(i+1)+". elembõl lett: "+row[i]);
			 }
			 if(row[i]==-0.0f) {
				 row[i] = 0;
			 }
			 
		 }
		 tempArray[rowIndex] = row;
		 System.out.println("");
	}
	
	/**
	 * Kiszámolja a b vektort.
	 * @param rowIndex
	 * @param columnIndex
	 */
	private void calcNewB(int rowIndex, int columnIndex) {
		float[] row = tempB;
		 for(int i = 0; i < row.length; i++) {
			 if (i == rowIndex) {
				 row[i] = row[i] * newGenValue;
				 System.out.println("Módosult B (gen.elem sorában) elem értéke: "+(i+1)+". elembõl lett: "+row[i]);
			 }
			 else{
				 row[i] = row[i]-(-b[rowIndex]*tempArray[i][columnIndex]);
				 System.out.println("Módosult B elemek értéke: "+(i+1)+". elembõl lett: "+row[i]);
			 }
			 if(row[i]==-0.0f) {
				 row[i] = 0;
			 }
		 }
		 tempB = row;
		 System.out.println("");
	}

	/**
	 * Kiszámolja a c vektort.
	 * @param rowIndex
	 * @param columnIndex
	 */
	private void calcNewC(int rowIndex, int columnIndex) {
		float[] col = tempC;
		 for(int i = 0; i < col.length; i++) {
			 if (i == columnIndex) {
				 col[i] = col[i] * (-1 * newGenValue);
				 System.out.println("Módosult C (gen.elem oszlopában) elem értéke: "+(i+1)+". elembõl lett: "+col[i]);
			 }
			 else{
				 col[i] = col[i]-(c[columnIndex]*tempArray[rowIndex][i]);
				 System.out.println("Módosult C elemek értéke: "+(i+1)+". elembõl lett: "+col[i]);
			 }
			 if(col[i]==-0.0f) {
					col[i] = 0;
			 }
		 }
		 tempC = col;
		 System.out.println("");
	}
	
	/**
	 * Kiszámolja a célfüggvény értékét
	 * @param rowIndex
	 * @param colIndex
	 */
	private void calcNewResult(int rowIndex, int colIndex) {
		float[] col = tempB;
		float[] row = c;
		result = result - row[colIndex] * col[rowIndex];
		System.out.println("Célfüggvény új értéke: "+result);
		System.out.println("");
	}
		 
	/**
	 * Kiszámolja a mátrix többi elemét.
	 * @param rowIndex
	 * @param colIndex
	 */
	private void calcNewMatrix(int rowIndex, int colIndex) {
		float[][] oldMatrix = A;
		float[][] newMatrix = tempArray;
		for (int i=0; i < newMatrix.length; i++){
			for (int j=0; j < newMatrix[i].length; j++){
				// a korábban kitöltött (generáló elem sora, oszlopa) mezõk figyelmen kívûl hagyása
				if (i == rowIndex || j == colIndex){
					continue;
				}
				// új koordináta = régi koordináta - (generáló elem sorának új értéke * generáló elem régi oszlopából a megfelelõ koordináta)
				newMatrix[i][j] = oldMatrix[i][j] - (newMatrix[rowIndex][j] * oldMatrix[i][colIndex]);
				if(newMatrix[i][j]==-0.0f) {
					newMatrix[i][j] = 0;
				}
			}
		}
		// mentés
		tempArray = newMatrix;
		System.out.println("");
	}

	/**
	 * Fájlból beolvasás, a mintaként kapott input alapján dinamikusan (n * m)
	 */
	private void loadFile() {
		FileInputStream fstream = null;
		BufferedReader br = null;
		try{
			// filePath - alapvetõen a fent megadott elérési út, futtatási paraméter megadása esetén pedig az
            fstream = new FileInputStream(filePath );
            br = new BufferedReader(new InputStreamReader(fstream));
            String strLine;
            N = Integer.parseInt(br.readLine());
            M = Integer.parseInt(br.readLine());
            br.readLine();
            
            // megfelelõ méretû tömbök létrehozása a tároláshoz
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
            
            // eredeti mátrix, vektorok elmentése (deep copy)
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
	 * Megfelelõ méretû tömbök létrehozása a tároláshoz.
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
	 * Kétdimenziós tömb teljes másolása egy újba.
	 * @param original - az eredeti 2d tömb
	 * @return a másolat
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
	 * Az A mátrix kiírása a konzolra.
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
	 * Paraméterben kapott vektor kiírása a konzolra.
	 * @param a - egész számokból álló vektor (u, v)
	 */
	private void printArray(int[] a) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i]+" ");
		}
		System.out.println("");
	}

	/**
	 * Paraméterben kapott vektor kiírása a konzolra.
	 * @param a - lebegõpontos számokból álló vektor (pl. b, c)
	 */
	private void printArray(float a[]) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i]+" ");
		}
		System.out.println("");
	}

	/**
	 * Visszaadja a pozitív c elemek listáját.
	 * Rendezi a c vektort, majd a választható értékeket beleírjük egy tömbbe.
	 * @return - index tömb
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
					System.out.println("A "+tempC[i]+" < 0, ezért biztosan nem választható a c - bõl");
				}
			} else {
				if(tempC[i] <= 0) {
					needDelete.add(i);
					System.out.println("A "+tempC[i]+" <= 0, ezért biztosan nem választható a c - bõl");
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
	 * Szûk keresztmetszet számítása 
	 * @param matrixIndex
	 * @return - legpozitívabb választás indexe
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
	 * Készí egy másolatot a kapott 2dimenziós tömb kívánt oszlopáról.
	 * @param matrixIndex - a kért oszlop indexe
	 * @param source - 2 dimenziós tömb
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
	 * Visszaadja a 2 dimenziós tömb kívánt sorát.
	 * @param source - 2 dimenziós tömb
	 * @param matrixIndex - a kért sor indexe
	 * @return
	 */
	private float[] getRow(float[][] source, int matrixIndex) {
		return source[matrixIndex];
	}

	/**
	 * Szûk keresztmetszethez szükséges osztás
	 * @param value
	 * @param index
	 * @return
	 */
	public float getDivisonResultForRow(float value, int index) {
		return b[index]/value;
	}
}
