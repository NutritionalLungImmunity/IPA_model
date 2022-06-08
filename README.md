# IPA model

To install, download the model using either `git` or by downloading and extracting a zip file. Then open a terminal and enter the model's directory. Run  
```
pipenv install
```
to install the model's dependencies in a virtual environment. (The only dependency is `numpy` at the time of writing.)  The virtual environment can be entered with
```
pipenv shell
```

To run, enter the virtual environment as above, then on linux or mac:
```
PYTHONPATH=. python edu/uf/main/Model.py
```
or on windows:
```
set PYTHONPATH=.
python edu/uf/main/Model.py
```

## Contact Information

<a href="mailto:henrique.deassis@medicine.ufl.edu">Henrique de Assis Ribeiro</a>, Laboratory for Systems Medicine, College of Medicine, University of Florida 

<a href="mailto:reinhard.laubenbacher@medicine.ufl.edu">Reinhard Laubenbacher</a>, Laboratory for Systems Medicine, College of Medicine, University of Florida