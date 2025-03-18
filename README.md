# SMILES.bot

SMILES.bot is a discord bot that renders Simplified Molecular Input Line Entry System strings into a image using `RDKit`.

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
