# SMILES.bot

SMILES.bot is a discord bot that renders Simplified Molecular Input Line Entry System strings into a image using `RDKit`.

## Functionality
### The `/render` Command
SMILES.bot is able to render molecules as well as chemical equations, using `/render mlcl` and `/render rxn` respectively. 

#### ID conversion
The `/render mlcl` function has a built in ID converter, meaning it can render a compound given in:
- IUPAC Names
- Common Names
- CAS Numbers
- InChIKey

<sub> if you get an error code, that compound may not be in the `CIRpy` library <sub>

### The `/set` Command

#### Prefix
You can set a custom prefix in your server using `/set prefix` allowing you to use both message commands as well as Discord slash commands.

#### Color
SMILES.bot has a custom palette designed to look good with Discord's ui, and if you disagree you can change the rener color of any element using `/set color <element> <color>`

##### Default Values
![smiles_bot_default](https://github.com/user-attachments/assets/1b95ddfd-64d2-491c-8f7a-ef2829bbdf66)

#### Background Color
The render background is automatically set to Discord's default darkmode but can be changed to whatever color using `/set bgcolor`
<p align="center">
<img src = "https://github.com/user-attachments/assets/4ac957c5-f07f-465e-89a4-a1f502258441" />
</p>

## SMILES
### SMILES Syntax
#### Bonds
` . ` Disconnected structures such as ionic bonds and multiple compounds in an image.\
` - `Single bonds (optional/usually omitted).\
` = ` Double bonds.\
` # ` Tripple bonds.\
` $ ` Quadruple bonds.\
` : ` Aromatic 'one and a half' bonds.\
` / ` Single bond (directional) use for cis-trans stereochemistry.\
` \ ` Single bond (directional) use for cis-trans stereochemistry.
#### Branching
`( ) `are used to denote branches, for example Isopropyl Alcohol is `CC(O)C`. Bond type is put inside the parentheses like `CCC(=O)O`.
#### Stereochemistry
` / and \ `are used for cis-trans stereochemistry, for example (E)-1,2-Dichloroethene is `Cl/C=C/Cl`. \
` @ and @@ ` are used to denote S (`@`) and R (`@@`) stereocenters.
#### Isotopes and charges
` [ and ] ` are used for charges and isotopes, for example Carbon-14 is [14C] and Sodium Chloride is `[Na+].[Cl-]`
## Functions
`/help` \
Show this help menu.

`/smileshelp`\
Quick guide to SMILES.

`/smiles <SMILES string>`\
Generate a molecular structure image from a SMILES string.

`/auto_detect <enable|disable>`\
Enable or disable auto-detection SMILES strings within `&...&`.

# DEV
create env: python -m venv .env
activate env: .\.env\Scripts\activate
