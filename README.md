<!-- PROJECT LOGO -->
<br />


<h3 align="center">Resonance Absorption Effects of Fertile Isotopes in Fusion Blankets</h3>



## About The Project

To-do

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->

## Getting Started

### Downloading our code
1. Download (clone) the repository:

   ```bash
   git clone https://github.com/patrickpark910/dt-fusion-illicit
   cd dt-fusion-illicit
   ```

2. Create and activate a virtual environment. If you have Windows, you probably should be running OpenMC through WSL:

   ```bash
   wsl # if on Windows
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install -U pip
   ```

3. Install your OpenMC into our virtual environment:

   ```bash
   pip install -e /path/to/your/openmc/folder
   ```

4. Install our Python dependencies:

   ```bash
   pip install -r requirements.txt
   ```



### Running our code

Once you have this set up, you can now run our code.  All our calculations are controlled via `main.py`. You can configure the run type, materials, and execution behavior using command-line flags:
   ```bash
   python main.py -r [run_type] -b [blankets] -i [isotopes] [flags]
   ```
Our code has three main functionalities: plot OpenMC geometries (`-r plot`), draw Miller shapes of a tokamak and compute volumes of each cell (`-r volume`), or calculate various tallies within the model (`-r tallies`). The full list of options is shown below:


   | **Flag** | **Full** | **Description**                                              | **Default** |
   | -------- | ------------ | ------------------------------------------------------------ | ----------- |
   | `-r`     | `--run_type` | Choose **one**: `tallies`, `volume`, or `plot`.              | `tallies`   |
   | `-b`     | `--blankets` | Specify one or more blankets (space-separated), e.g., `-b FLiBe` or  `-b FLiBe DCLL HCPB` | All         |
   | `-i`     | `--isotopes` | Specify one or more fertile isotopes (space-separated), e.g., `-i U238` | All         |
   |          | `--no_xml`   | Runs through the code but doesn't print `model.xml` (useful for debugging) |             |
   |          | `--no_run`   | Doesn't execute the OpenMC calculation                       |             |
   |          | `--no_debug` | Turns off the debug print statements                         |             |

NB. The engine, by default, automatically prints XMLs, runs OpenMC, and prints debug statements. Use the flags `--no_xml --no_run --no_debug` to turn off the respective function.



### Usage Examples

Run tally calculations for a FLiBe blankets for all fertile isotopes (U238, Th232):
```bash
python main.py -r tallies -b FLiBe
```
Run tally calculations for FLiBe, DCLL blankets and only add U238:
```bash
python main.py -r tallies -b FLiBe DCLL -i U238 Th232
```
Run tally simulations for a HCPB blanket and only add Th232. Say I want to test what the volume ratios of Th232 to Li4SiO4 are, for which I have written some debug print statements, and I don't need to print the XML nor run OpenMC. 
```bash
python main.py -r tallies -b HCPB -i Th232 --no_xml --no_run
```
Generate the XML for a geometry plots without actually running the OpenMC plotting calculation:
```bash
python main.py -r plot --no_run
```


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue or shoot us an email. 

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->

## Contact

Emma R. Zoccoli - ezoccoli3@gatech.edu
Greta X. Li - gl3210@princeton.edu 
Patrick J. Park - pjp2136@columbia.edu 
Robert J. Goldston - Principal Investigator - rjg@princeton.edu 

Project Link: [https://github.com/patrickpark910/dt-fusion-illicit](https://github.com/patrickpark910/dt-fusion-illicit)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->

## Acknowledgments

* []()
* []()
* []()

<p align="right">(<a href="#readme-top">back to top</a>)</p>
