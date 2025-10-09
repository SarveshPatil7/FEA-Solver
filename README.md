# FEA-Solver

Code cleanup in progress.

This is a hands on project for thermal FEA for 2D files. Only supports DC2D4 elements at the moment.

---

## Comparison of results with Abaqus

Testing the performance with distorted elements
<p float="left">
  <figure>
    <img src="Images/distort_abaq.png" width="380"/>
    <figcaption align="center"><b>Results from Abaqus</b></figcaption>
  </figure>
  <figure>
    <img src="Images/distort_py.png" width="380"/>
    <figcaption align="center"><b>Results from the custom solver in python</b></figcaption>
  </figure>
</p>

Results from a higher element count simulation
<p float="left">
  <figure>
    <img src="Images/notched_abaq.png" width="380"/>
    <figcaption align="center"><b>Results from Abaqus</b></figcaption>
  </figure>
  <figure>
    <img src="Images/notched_py.png" width="380"/>
    <figcaption align="center"><b>Results from the custom solver in python</b></figcaption>
  </figure>
</p>

---

## Usage

