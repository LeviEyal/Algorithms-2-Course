<div dir='rtl' lang='he'> 

<h1 align="center">פסאודו קודים - קורס אלגוריתמים 2 סמסטר ב' 2021</h1>

- ## [תוכן עניינים](#------------)
    - [פלויד וורשאל - Floyd Warshall](#פלויד-וורשאל-Floyd-Warshall)
    - [תת מערך עם סכום תאים מקסימלי](#תת-מערך-עם-סכום-תאים-מקסימלי)
    - [תת מטריצה עם סכום תאים מקסימלי](#תת-מטריצה-עם-סכום-תאים-מקסימלי)
    - [עץ פורש מינימלי](#עץ-פורש-מינימלי)
        - [קרוסקל](#עץ-פורש-מינימלי)
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

# פלויד וורשאל - Floyd Warshall
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
    new empty-string paths[N,N]
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
</div>
	
```python
Has-Negative-Cycles-UnDirected(g[N,N]) : O(n^2)
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
Degrees-of-vertices(g[N,N]) : O(n^2)
    new D[N] = {0,0,..,0}
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
Degrees-of-vertices(g[N,N]) : O(n^2)
    new D[N] = {0,0,..,0}
    for i=0 to N :
        for j=0 to N :
            if g[i,j] = true:
                D[i]++
                D[j]++
    sort(D)
    return D
```
<div dir='rtl' lang='he'>

# תת מערך עם סכום תאים מקסימלי	

* גרסה מקוצרת כשצריך רק את הסכום ללא אינדקסים
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

* גרסה רגילה כשצריך את הסכום ואת האינדקסים של תת המערך
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

גרסה מעגלית:
* סיבוכיות: `O(n)`
</div>

```python
Best-Cycle(A[N]) :
    totalSum = sum(A)
    return max(Best(A), totalSum - Best(-A))
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
        new help[cols]
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
        new help[cols]
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
        new node z
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
        new node z
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
    for each v in V :
        makeSet(v)
    sort E by weights
    for each (u,v) in E :
        if findSet(u) ≠ findSet(v) :
            T.add((u,v))
            union(u, v)
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
    for each v in V :
        makeSet(v)

    while |T| < |V|-1 :
        new cheapest[|V|] := array of edges
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
    for each v in V :
        v.color = WHITE
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
	
מציאת מספר רכיבי קשירות בגרף:
* סיבוכיות: `O(|V|+|E|)`
</div>

```python
DFS-Number-Of-Connected-Components(G) :
    counter = 0
    for each v in V :
        v.color = WHITE
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
    for each v in V :
        v.color = WHITE
    for each v in V :
        if v.color = WHITE :
            DFS-VISIT(G, v)
    return false
```

	
```python
DFS-VISIT(G, n) :
    n.color = GRAY
    for each v in adj[n] :
        if v.color = GRAY : 
            return true
        if v.color = WHITE :
            DFS-VISIT(G, v)
    n.color = BLACK
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
            C.add(v)
            stack.pop()
        else :
            u = first of adj[v]
            stack.push(u)
            G.remove(u,v)
            G.remove(v,u)
            u.degree--
            v.degree--
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
            for each v in adj[leaf] :
                if --v.deg = 1 :
                    future.add(v)
                n--
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
* סיבוכיות: `O((V+E)logV)`
	
</div>

```python
Dijkstra(G, s) :
    for v in V :
        v.dist = ∞

    s.dist = 0
    Q = V
    while Q ≠ ∅ :
        u = Q.extractMin() //by dists
        for v in adj[u] :
            if v.dist > u.dist + w(u,v) :
                v.dist = u.dist + w(u,v)
                v.prev = u

```


