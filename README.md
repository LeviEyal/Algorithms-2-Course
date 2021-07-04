<div dir='rtl' lang='he'> 

<h1 align="center">פסאודו קודים - קורס אלגוריתמים 2 סמסטר ב' 2021</h1>

- ## [תוכן עניינים](#------------)
    - [פלויד וורשאל - Floyd Warshall](#Floyd-Warshall)
    - [בעיית הבקבוקים](#בעיית-הבקבוקים)
    - [תת מערך עם סכום תאים מקסימלי](#תת-מערך-עם-סכום-תאים-מקסימלי)
        - [בעיית תחנות הדלק](#תת-מערך-עם-סכום-תאים-מקסימלי)
    - [תת מטריצה עם סכום תאים מקסימלי](#תת-מטריצה-עם-סכום-תאים-מקסימלי)
    - [עץ פורש מינימלי](#עץ-פורש-מינימלי)
        - [קרוסקל](#עץ-פורש-מינימלי)
        - [קרוסקל הפוך](#עץ-פורש-מינימלי)
        - [פרים](#עץ-פורש-מינימלי)
        - [בורובקה](#עץ-פורש-מינימלי)
    - [מעבר על גרפים](#מעבר-על-גרפים)
        - [BFS](#מעבר-על-גרפים)
        - [DFS](#מעבר-על-גרפים)
        - [מציאת מספר רכיבי קשירות בגרף](#מעבר-על-גרפים)
        - [בדיקת קיום מעגל בגרף](#מעבר-על-גרפים)
        - [בדיקת קיום מעגל אויילר בגרף](#מעבר-על-גרפים)
        - [מציאת מעגל אויילר בגרף](#מעבר-על-גרפים)
    - [קוטר ורדיוס בעצים](#קוטר-ורדיוס-בעצים)
        - [מציאת קוטר, רדיוס ומרכזים של עץ ע"י אלגוריתם שריפת עלים](#קוטר-ורדיוס-בעצים)
        - [מציאת קוטר של עץ ע"י אלגוריתם סריקה לעומק](#קוטר-ורדיוס-בעצים)
        - [מציאת המסלול של קוטר בעץ ע"י אלגוריתם סריקה לעומק
](#קוטר-ורדיוס-בעצים)
        - [בדיקה האם קודקוד נתון נמצא על קוטר בעץ
](#קוטר-ורדיוס-בעצים)
    - [קוד האפמן - Huffman code](#קוד-האפמן)
        - [בניית עץ](#קוד-האפמן)
        - [הוצאת הייצוגים הבינאריים של כל תו מהעץ](#קוד-האפמן)
        - [קוד האפמן - כאשר מערך התוים כבר ממוין](#קוד-האפמן)
    - [דייקסטרה Dijkstra](#דייקסטרה)
    - [בניית עץ מרשימת דרגות](#בניית-עץ-מרשימת-דרגות)
    - [איזומורפיזם של עצים](#איזומורפיזם-של-עצים)


# Floyd Warshall
גרף לא ממושקל: מפעילים כדי לבדוק אם קיים מסלול כלשהו בין קודקוד מסוייים לקודקוד אחר 
* סיבוכיות: `O(n^3)`
</div>

```python
FW-UnWheighted(g[N,N]) :
    for k=0 to N :
        for i=0 to N :
            for j=0 to N :
           		g[i,j] = g[i,j] || (g[i,k] AND g[k,j])
```
<div dir='rtl' lang='he'>
    
גרף ממושקל - מציאת כל המרחקים הקצרים ביותר בגרף בין כל זוג קודקודים:
* סיבוכיות: `O(n^3)`
</div>

```python
FW-wheighted(g[N,N]) :
    for k=0 to N :
        for i=0 to N :
            for j=0 to N :
                if  g[i,j] > g[i,k] + g[k,j] :
                    g[i,j] = g[i,k] + g[k,j]
```
<div dir='rtl' lang='he'>
גרף ממושקל - מציאת כל המסלולים הקצרים ביותר בגרף בין כל זוג קודקודים מיוצגים ע"י מחרוזות:

* סיבוכיות: `O(n^3)`
</div>

```python
FW-wheighted-paths(g[N,N]) :
    New empty-string paths[N,N]
    for i=0 to N :
        for j=0 to N :
            paths[i,j] = i + "->" + j

    for k=0 to N :
        for i=0 to N :
            for j=0 to N :
                if  g[i,j] > g[i,k] + g[k,j] :
                    g[i,j] = g[i,k] + g[k,j]
                    paths[i,j] = paths[i,k] + "->" + paths[k,j]
    return paths
```

<div dir='rtl' lang='he'>

גרף לא מכוון - בדיקת קשירות:
* מספיק להפעיל פלוייד וורשל על הגרף ואז לבדוק את השורה הראשונה. כי אם מהקודקוד הראשון יש מסלול לכולם אז אפשר ישר להגיד שהגרף קשיר 
* סיבוכיות: `O(n^3)`
</div>

```python
isConnected-UnDirected-graphs(g[N,N]) :
    FW(g)
    for i=1 to N :
        if g[0,i] = false :
            return false
    return true

```
<div dir='rtl' lang='he'>

בדיקת קשירות של גרפים לא מכוונים:
* סיבוכיות: `O(n^3)`
</div>
	
```python
isConnected-Directed-graphs(g[N,N]) :
    FW(g)
    for i=0 to N :
        if g[0,i] = ∞ :
            return false
    return true
```
<div dir='rtl' lang='he'>

בדיקה האם קיים מעגל שלילי בגרף מכוון:

* מספיק להסתכל על האלכסון הראשי אחרי הפעלת `Floyd-Washall` 
* סיבוכיות: `O(n^3)`
</div>
	
```python
Has-Negative-Cycles-Directed(g[N,N]) :
    FW(g)
    for i=0 to N :
        if g[i,i] < 0 :
            return true
    return false
```
<div dir='rtl' lang='he'>
בדיקה האם קיים מעגל שלילי בגרף לא מכוון:

* אין צורך להפעיל `Floyd-Washall`
* אם קיימת צלע עם משקל שלילי אז קיים מעגל שלילי בגרף לא מכוון
* סיבוכיות: `O(n^2)`
</div>
	
```python
Has-Negative-Cycles-UnDirected(g[N,N]) : 
    for i=0 to N :
        for j=0 to N :
            if g[i,j] < 0:
                return true
    return false
```
<div dir='rtl' lang='he'>

חישוב דרגות קודקודים ממטריצת מרחקים בגרף ממושקל ולא מכוון:
* סיבוכיות: `O(n^2)`
</div>

```python
Degrees-of-vertices(g[N,N]) :
    New D[N] = {0,0,..,0}
    for i=0 to N :
        for j=0 to N :
            if i ≠ j AND g[i,j] ≠ ∞:
                D[i]++
                D[j]++
    sort(D)
    return D
```
<div dir='rtl' lang='he'>

חישוב דרגות קודקודים ממטריצת שכנויות בגרף לא ממושקל ולא מכוון:
* סיבוכיות: `O(n^2)`
</div>

```python
Degrees-of-vertices(g[N,N]) :
    New D[N] = {0,0,..,0}
    for i=0 to N :
        for j=0 to N :
            if g[i,j] = true:
                D[i]++
                D[j]++
    sort(D)
    return D
```
<div dir='rtl' lang='he'>

# בעיית הבקבוקים
* סיבוכיות: `O(n^3)`
</div>

```python
Bottles-Problem(b1, b2) :
    mat = Bottles-Problem-Create-Matrix(b1, b2)
    return Bottles-Problem-Get-Paths(mat, max(b1, b2))
```
<div dir='rtl' lang='he'>

שלב ראשון -  בניית מטריצת המעברים ממצב למצב:
* סיבוכיות: `O(n^2)`
</div>

```python
Bottles-Problem-Create-Matrix(b1, b2) :
    size = (b1+1)*(b2+1)
    max = max(b1, b2)
    new mat[size, size]

    for i=0 to b1 :
        for j=0 to b2 :
            index = getindex(max, i, j)

            1. mat[index, getIndex(max, 0, j)] = 1
            2. mat[index, getIndex(max, i, 0)] = 1
            3. mat[index, getIndex(max, b1, j)] = 1
            4. mat[index, getIndex(max, i, b2)] = 1
            5. mat[index, getIndex(max, min(b1, i+j), i+j-min(b1, i+j))] = 1
            6. mat[index, getIndex(max, i+j-min(b2, i+j) , min(b2, i+j))] = 1
    return mat

getIndex(n, i, j) :
    return (n+1) * i + j
```
<div dir='rtl' lang='he'>

שלב שני - יצירת מטריצת מסלולים ממצב למצב:
* סיבוכיות: `O(n^2)`
</div>

```python
Bottles-Problem-Get-Paths(mat[size, size], n) :
    New paths[size, size]
    for i=0 to size :
        x1, y1 = getPairFromIndex(n, i)
        for j=0 to size :
            x2, y2 = getPairFromIndex(n, i)
            if mat[i, j] ≠ ∞ :
                paths[i, j] = "("+x1+","+y1+") → ("+x2+","+y2+")"
    for k=0 to size :
        for i=0 to size :
            for j=0 to size :
                if mat[i,j] > mat[i,k] + mat[k,j] :
                   mat[i,j] = mat[i,k] + mat[k,j]
                   paths[i,j] = paths[i,k] + "→" + paths[k,j]
    return paths
    
getPairFromIndex(n, i) :
    return { i/(n+1), i%(n+1) }
```

<div dir='rtl' lang='he'>

# תת מערך עם סכום תאים מקסימלי	

גרסה מקוצרת כשצריך רק את הסכום ללא אינדקסים
* סיבוכיות: `O(n)`
</div>
	
```python
Best-basic(A[N]) :
    maxSum = sum = -∞

    for i=0 to N :
        sum += A[i]
        maxSum = max(maxSum, sum)
        sum = max(max, 0)
    return max
```
<div dir='rtl' lang='he'>

גרסה רגילה כשצריך את הסכום ואת האינדקסים של תת המערך
* סיבוכיות: `O(n)`
</div>

```python
Best(A[N]) :
    max = sum = -∞
    start = end = 0

    for i=0 to N :
        sum += A[i]
        if max < sum :
            max = sum
            end = i
        if max < 0:
            sum = 0
            start = i+1
    return {max, start, end}
```
<div dir='rtl' lang='he'>

גרסה מעגלית כשלא צריך אינדקסים:
* סיבוכיות: `O(n)`
</div>

```python
Best-Cycle(A[N]) :
    totalSum = sum(A)
    return max(Best(A), totalSum - Best(-A))
``` 

<div dir='rtl' lang='he'>

גרסה מעגלית כולל אינדקסים:
* סיבוכיות: `O(n)`
</div>

```python
Best-Cycle(A[N]) :
    totalSum = sum(A)
    max1, start1, end1 = Best(A)
    max2, start2, end2 = Best(-A)
    
    if max1 > totalSum + max2 :
        return {max1, start1, end1}
    else
        return {totalSum + max2, end2 + 1, start2 - 1}
``` 

<div dir='rtl' lang='he'>

בעיית תחנות הדלק:
* מטרה למצוא תחנה ממנה אפשר ליסוע ולעשות סיבוב שלם.
* המערך A מייצג את כמות הדלק בכל תחנהת המערך B מייצג את עלות הנסיעה בדלק מתחנה לתחנה.
* סיבוכיות: `O(n)`
</div>

```python
Gas-Stations-Problem(A[N], B[N]) :
    C[N]
    sum = 0
    for i=0 to N :
        C[i] = A[i] - B[i]
        sum += C[i]
    
    if sum < 0 :
        return "No solution"
    return Best-Cycle(C)
``` 

<div dir='rtl' lang='he'>
	
# תת מטריצה עם סכום תאים מקסימלי
	
גרסה מקוצרת כשצריך רק את הסכום ללא אינדקסים:
* סיבוכיות: `O(n^3)`
</div>

```python
Super-Best-basic(A[rows,cols]):
    maxSum = 0
    for i=0 to rows :
        New help[cols]
        for j=i to rows :
            for k=0 to cols :
                help[k] += A[j,k]
            maxSum = max(maxSum, Best-basic(help))
    return maxSum

```
<div dir='rtl' lang='he'>

גרסה רגילה כשצריך את הסכום ואת האינדקסים של תת המטריצה:
* סיבוכיות: `O(n^3)`
</div>

```python
Super-Best(A[rows,cols]) :
    maxSum, start_x, start_y ,end_x, end_y = 0

    for i=0 to rows :
        New help[cols]
        for j=i to rows :
            for k=0 to cols :
                help[k] += A[j,k]
            tempSum, tempStart, tempEnd = best(help)
            if maxSum < tempSum :
                maxSum = tempSum
                start_x = tempStart
                end_x = tempEnd
                start_y = i
                end_y = j
    return {maxSum, start_x, start_y, end_x, end_y}

```
<div dir='rtl' lang='he'>

# קוד האפמן

שלב ראשון בניית העץ:
* סיבוכיות: `O(n*log(n))`	
</div>

```python
Huffman(C) :
    Q = C
    while |Q| > 1 :
        New node z
        z.left = x = q.extractMin()
        z.right = y = q.extractMin()
        z.freq = x.freq + y.freq
        Q.insert(z)
    return Q.extrctMin()
```
<div dir='rtl' lang='he'>

שלב שני - הוצאת הייצוגים הבינאריים של כל תו מהעץ:
* סיבוכיות: `O(n)`	
</div>
	
```python
Huffman-process(root = Huffman(C)) :
    if(root.left == null AND root.right == null) :
        print(root.char)
    else:
        print('0' + Huffman-process(root.left)
        print('1' + Huffman-process(root.right)

```

<div dir='rtl' lang='he'>
	
קוד האפמן - כאשר מערך התוים ממוין לפי שכיחויות:
* סיבוכיות: `O(n)`	
</div>
	
```python
Huffman-sorted(C) :
    Q1 = C
    Q2 = ∅
    while |Q1| + |Q2| > 1 :
        New node z
        z.left = x = getMin(Q1,Q2)
        z.right = y = getMin(Q1,Q2)
        z.freq = x.freq + y.freq
        Q2.insert(z)
    return getMin(Q1, Q2)

getMin(Q1,Q2):
    if Q1 = ∅ :
        return Q2.dequeue()
    if Q2 = ∅ :
        return Q1.dequeue()
    if Q1.front < Q2.front :
        return Q1.dequeue()
    else return Q2.dequeue()
```
<div dir='rtl' lang='he'>
	
# עץ פורש מינימלי

יצירת עץ פורש מינימלי - קרוסקל:
* סיבוכיות: `O(|E|log|E|)`
	
</div>
	
```python
MST-Kruskal(G) :
    T = ∅
    for each v in V : makeSet(v)
    sort E by weights
    for each (u,v) in E :
        if findSet(u) ≠ findSet(v) :
            T.add((u,v))
            union(u, v)
    return T

```
<div dir='rtl' lang='he'>
	
יצירת עץ פורש מינימלי - קרוסקל הפוך:
* הרעיון הוא להתחיל עם "עץ" שיש בו את כל הצלעות של G. נעבור על הצלעות מהגדולות לקטנות ונשאל על כל צלע האם אפשר לנתק אותה מהגרף והוא עדיין יישאר קשיר. אם כן נמחק את הצלע מהעץ, ככה באותו האופן עד שנגיע לעץ בגודל תקין.
* סיבוכיות: `O(|E|(|E|+|V|))`
	
</div>
	
```python
MST-Reversed-Kruskal(G) :
    T = E
    Q = E
    while |T| > |V|-1 :
        e = Q.extractMax()
        G.removeEdge(e)
        if isConnected = false :
            G.addEdge(e)
    return T
```
<div dir='rtl' lang='he'>
	
יצירת עץ פורש מינימלי - פרים:
* סיבוכיות: `O(|V|+|E|log|V|)`
	
</div>

```python
MST-Prim(G) :
    T = ∅
    for each v in V :
        v.key = ∞
        v.prev = null
    V[0].key = 0
    Q = V
    while Q ≠ ∅ :
        u = Q.extractMin()
        if u.prev ≠ null :
            T.add((u, u.prev))
        for each v in adj[u] :
            if v in Q AND v.key > w(u,v) :
                v.key = w(u,v)
                v.prev = u
    return T

```
<div dir='rtl' lang='he'>
	
יצירת עץ פורש מינימלי - בורובקה:
* סיבוכיות: `O(|E|log|V|)`
	
</div>

```python
MST-Boruvka(G) :
    T = ∅
    for each v in V : makeSet(v)

    while |T| < |V|-1 :
        New cheapest[|V|] := array of edges
        for (u,v) in E :
            g1 = findSet(u)
            g2 = findSet(v)
            if g1 ≠ g2 :
                if w(cheapest[g1]) > w(u,v) :
                    w(cheapest[g1]) = w(u,v)
                if w(cheapest[g2]) > w(u,v) :
                    w(cheapest[g2]) = w(u,v)
        for i=0 to |v| :
            if cheapest[i] ≠ null :
                T.add(cheapest[i])
                union(cheapest[i].u, cheapest[i].v)
    return T

```
<div dir='rtl' lang='he'>
	
# מעבר על גרפים
## BFS
מעבר על גרף באמצעות סריקה לרוחב:
* סיבוכיות: `O(|V|+|E|)`
</div>
	
```python
BFS(G, s) :
    for each v in V :
        v.color = WHITE, v.dist = ∞, v.prev = null

    s.color = GRAY, s.dist = 0, s.prev = -1
    Q = {s}
    while Q ≠ ∅ :
        u = Q.dequeue()
        for each v in adj[u] :
            if v.color = WHITE :
                v.color = GRAY
                v.dist = u.dist + 1
                v.prev = u
                Q.enqueue(v)
        u.color = BLACK
```

<div dir='rtl' lang='he'>

## DFS
מעבר על גרף באמצעות סריקה לעומק:
* סיבוכיות: `O(|V|+|E|)`
</div>
	
```python
DFS(G) :
    for each v in V : v.color = WHITE
    for each v in V :
        if v.color = WHITE :
            DFS-VISIT(G, v)

DFS-VISIT(G, n) :
    n.color = GRAY
    for each v in adj[n] :
        if v.color = WHITE :
            DFS-VISIT(G, v)
    n.color = BLACK
```
<div dir='rtl' lang='he'>
	
בדיקה האם צלע נתונה נמצאת על מעגל:
* סיבוכיות: `O(|V|+|E|)`
</div>

```python
Is-Edge-On-Cycle(G, e) :
    u = e.u
    v = e.v
    G.removeEdge(e)
    BFS(v)
    return (u.color ≠ WHITE)
```
<div dir='rtl' lang='he'>
	
מציאת מספר רכיבי קשירות בגרף:
* סיבוכיות: `O(|V|+|E|)`
</div>

```python
DFS-Number-Of-Connected-Components(G) :
    counter = 0
    for each v in V : v.color = WHITE
    for each v in V :
        if v.color = WHITE :
            counter++
            DFS-VISIT(G, v)
    return counter
```
<div dir='rtl' lang='he'>
	
בדיקה האם קיים מעגל בגרף:
* סיבוכיות: `O(|V|+|E|)`
</div>
	
```python
DFS-Has-Cycle(G) :
    for each v in V : v.color = WHITE
    for each v in V :
        if v.color = WHITE :
            if DFS-VISIT(G, v) = true :
                return true
    return false
```

	
```python
DFS-VISIT(G, n) :
    n.color = GRAY
    for each v in adj[n] :
        if v.color = GRAY : 
            return true
        if v.color = WHITE :
            if DFS-VISIT(G, v) = true :
                return true
    n.color = BLACK
    return false
```
<div dir='rtl' lang='he'>
	
בדיקה האם קיים מעגל אויילר בגרף:
* סיבוכיות: `O(|V|+|E|)`
</div>
	
```python
Has-Euler-Cycle(G) :
    for each v in V :
        if v.degree % 2 ≠ 0 :
            return "No euler Cycle detected"
    if isConnected(G) :
        return Euler-Cycle(G)
```

<div dir='rtl' lang='he'>
	
מציאת מעגל אויילר בגרף:
* סיבוכיות: `O(|V|+|E|)`
</div>
	
```python
Euler-Cycle(G) :
    select some vertex s from G 
    list C = ∅
    stack = {s}
    while stack ≠ ∅ :
        v = stack.peek()
        if v.degree = 0 :
            C.add(stack.pop())
        else :
            u = first of adj[v]
            stack.push(u)
            G.remove(u,v)
            G.remove(v,u)
    return C        
```

<div dir='rtl' lang='he'>
	
בדיקת קשירות בגרף ע"י סריקה לרוחב:
* סיבוכיות: `O(|V|+|E|)`
</div>
	
```python
BFS-isConnected(G) :
    select some vertex s from G 
    BFS(G, s)
    for each v in V :
        if v.color ≠ BLACK :
            return false
    return true 
```

<div dir='rtl' lang='he'>
	
בדיקה האם קיים מסלול אויילר בגרף:
* סיבוכיות: `O(|V|+|E|)`
</div>

```python
Euler-HasEuler-path(G) :
    counter = 0
    for each v in V :
        if v.degree % 2 ≠ 0 :
            counter++
    return isConnected(G) = true AND counter = 2

```
<div dir='rtl' lang='he'>

# קוטר ורדיוס בעצים
מציאת קוטר, רדיוס ומרכזים של עץ ע"י אלגוריתם שריפת עלים:
* סיבוכיות: `O(|V|+|E|)`
	
</div>

```python
Diameter-Fire(T) :
    n = |V|, radius = 0, leaves = ∅
    for each v in V :
        if v.deg = 1 :
            leaves.add(v)
    
    while n > 2 :
        future = ∅
        for each leaf in leaves :
            leaf.deg = 0
            n--
            for each v in adj[leaf] :
                if --v.deg = 1 :
                    future.add(v)
        radius++
        leaves = future

    diameter = radius*2 + |leaves|-1
    centers[] = leaves
    return {centers, radius, diameter}
```

<div dir='rtl' lang='he'>

מציאת קוטר של עץ ע"י אלגוריתם סריקה לעומק:
* סיבוכיות: `O(|V|+|E|)`
	
</div>

```python
Diameter-DFS(T) :
    maxDist = 0, MaxDistVertex = 0 
    DFS(T, 0)
    for each v in V:
        if maxDist < v.dist :
            maxDist = v.dist
            MaxDistVertex = v 

    DFS(T, MaxDistVertex)

    for each v in V:
        if maxDist < v.dist :
            maxDist = v.dist
            
    return maxDist
```

<div dir='rtl' lang='he'>

מציאת המסלול של קוטר בעץ ע"י אלגוריתם סריקה לעומק DFS:
* סיבוכיות: `O(|V|+|E|)`
	
</div>

```python
Diameter-Path-DFS(T) :
    maxDist = 0, MaxDistVertex = 0 
    DFS(T, 0)
    for each v in V:
        if maxDist < v.dist :
            maxDist = v.dist
            MaxDistVertex = v 

    DFS(T, MaxDistVertex)

    for each v in V:
        if maxDist < v.dist :
            maxDist = v.dist
            MaxDistVertex = v

    path = ∅
    while MaxDistVertex != -1 :
        path.add(MaxDistVertex)
        MaxDistVertex = MaxDistVertex.prev

    return path
```
<div dir='rtl' lang='he'>

בדיקה האם קודקוד נתון נמצא על קוטר בעץ:
* הבעיה דומה לבדיקה האם נקודה נמצאת על צלע מסוים.
* שלב ראשון: נמצא את הקוטר (ע"י אלגוריתם שריפה או באמצעות `DFS`) 
* שלב שני: נמצא את הקודקוד הכי רחוק מקודקוד `v` - מובטח שהוא קצה אחד של הקוטר. 
* שלב שלישי: נמצא את הצלע הראשונה מקודקוד `v` במסלול לעבר הקודקוד שמצאנו בשלב הקודם. נשמור את המרחק.
* שלב רביעי: נסיר את הצלע הזאת מה שבעצם מסיר את כל הענף
* שלב חמישי: נעשה שוב את שלב 2. מובטח שהפעם הקודקוד הכי רחוק מ `v` הוא הקצה השני של הקוטר.
* שלב אחרון: אם הקוטר שווה למרחק של `v` מהקצה הראשון של הקוטר ועוד מרחקו מהקצה השני - אזי הקודקוד נמצא על הקוטר (אחד מהם במידה ויש כמה)
* סיבוכיות: `O(|V|+|E|)`
	
</div>

```python
is-Vertex-on-Diameter(T, v) :
    diameter = Diameter(T)
    BFS(v)
    maxDist1 = maxDistIndex = 0
    for each v in V :
        if maxDist1 < v.dist :
            maxDist1 = v.dist
            maxDistIndex = v

    while maxDistIndex.prev ≠ v :
        maxDistIndex = maxDistIndex.prev

    remove((v, maxDistIndex))
    BFS(v)
    maxDist2 = 0
    for each v in V :
        if maxDist2 < v.dist :
            maxDist2 = v.dist
    
    return (diameter == maxDist1 + maxDist2)

```

<div dir='rtl' lang='he'>

# דייקסטרה
מציאת המרחקים הקצרים ביותר בגרף מקודקוד נתון:
* סיבוכיות: `O((V+E)logV)`  (כאשר משתמשים בערימה בינארית)
	
</div>

```python
Dijkstra(G, s) :
    for v in V : v.dist = ∞

    s.dist = 0
    Q = V
    while Q ≠ ∅ :
        u = Q.extractMin() //by dists
        for v in adj[u] :
            if v.dist > u.dist + w(u,v) :
                v.dist = u.dist + w(u,v)
                v.prev = u

```
<div dir='rtl' lang='he'>

# בניית עץ מרשימת דרגות
* צריך לבדוק קודם שאכן מתקיים  `Sum(degs) = 2*(|V|-1)`  לפני שמתחילים את האלגוריתם.
* סיבוכיות: `O(VlogV)` אם המערך לא ממוין. `O(V)` אם המערך ממוין כבר.
* הערה: האינדקס של האיבר הראשון במערך הוא 0, האינדקס של האיבר האחרון במערך הוא N-1.
	
</div>

```python
GenerateTreebyDegrees(deg[N]) :
    sort(deg)

    j = 0
    while deg[j] = 1 : j++
    
    New Tree[N,N] = {false}
    for i=0 to N-2 :
        Tree[i,j] = true
        Tree[j,i] = true
        if --deg[j] = 1 :
            j++
    Tree[N-2,N-1] = true
    Tree[N-1,N-2] = true
    return Tree
```

<div dir='rtl' lang='he'>

# איזומורפיזם של עצים
* הרעיון הוא ליצור ייצוג בינארי כמחרוזת לכל עץ. אם המחרוזת של העץ הראשון זהה למחרוזת של העץ השני, אזי שני העצים איזומורפיים.
* יצירת הייצוג הבינארי כמחרוזת לעץ מתחילה רקורסיבית מהשורש כלפי מטה, כאשר כל עלה מיוצג ע"י "01" וכל אבא משרשר את הייצוגים של הילדים שלו בסדר ממוין כאשר משמאל יש "0" ומימין יש "1".
* אם לא נתון מה השורש של העץ אז נשתמש באלגוריתם שריפה למציאת מרכז העץ והוא יהיה השורש.
* אם לעץ הראשון יש מרכז אחד ולשני יש שני מרכזים, אזי העצים לא איזומורפיים.
* אם יש 2 מרכזים לעץ, נבדוק אם המחרוזת של העץ הראשון שווה למחרוזת של העץ השני בחילופי תפקידים בין המרכזים שנמצאו.
* סיבוכיות: `O(VlogV)`
	
</div>

```python
AHU-Tree-Isomorphism(T1, T2) :
    r1 = T1.root
    r2 = T2.root
    if r1 = null AND r2 = null :
        r1 = findCenter(T1)    // by Fire Algorithm
        r2 = findCenter(T2)
    New global List<String> childrenCodes[|V|]
    code1 = findCode(r1)
    code2 = findCode(r2)
    return (code1 == code2)
```

```python
findCode(u)
    u.color = BLACK
    if u is a leaf :
        return "10"

    for each v in adj[u] :
        if v.color = WHITE :    // then v is child of u
            childrenCodes[u].add(findCode(v))
    Sort(childrenCodes[u])
    temp = ""
    for each s in childrenCodes[u] :
        temp += s
    return "1"+temp+"0"
```






