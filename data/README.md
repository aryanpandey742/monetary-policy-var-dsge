# Data Sources

Raw data is not included in this repository. To reproduce the analysis, you need an Excel
workbook with two sheets — `Quarterly` and `Monthly` — each with an `observation_date` column
plus the series below, all of which are freely available from the Federal Reserve Economic
Data (FRED) database at https://fred.stlouisfed.org/.

**Quarterly sheet** (1966 Q1 onward):

| Column | FRED series |
|---|---|
| `GDPC1` | Real Gross Domestic Product |
| `HOANBS` | Nonfarm Business Sector: Hours Worked for All Employed Persons |
| `FDEFX` | Federal Government: National Defense, Consumption Expenditures |
| `USAGDPDEFQISMEI` | GDP Deflator |
| `DFF` | Federal Funds Effective Rate |

**Monthly sheet**:

| Column | FRED series |
|---|---|
| `CLF16OV` | Civilian Labor Force Level |

Place the workbook in this folder and update the filename referenced in the notebook to match.
