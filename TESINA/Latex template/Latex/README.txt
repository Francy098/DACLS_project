Questa cartella contiene una serie di file che permette di scrivere la tesina/tesi in ambiente Latex.
Nel file tesi.cls sono contenute le definizioni dei vari ambienti.
Nel file tesi.tex è presente un esempio di utilizzo della classe tesi.cls.

Per introdurre i dati nel frontespizio è necessario settare le variabili author, corso, titolo, relatore, correlatoretxt, 
correlatore, annoaccademico come mostrato nell'esempio tesi.tex.

Se il correlatore non è presente sostituire i campi correlatore e correlatoretxt con i seguenti valori:
\correlatoretxt{}
\correlatore{\\} 

Nel caso non si volesse utilizzare la seguente classe si chiede comunque di rispettare le specifiche sul 
frontespizio, sul carattere e sui margini in essa presenti. Si richiede la scrittura dell'abstract.

Il file tesi.cls si basa sul file unthesis.cls scaricabile dal sito
http://www.ctan.org/tex-archive/help/Catalogue/entries/ucthesis.html#Download