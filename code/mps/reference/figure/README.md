# Initial-velocity MPS: revised six-panel figure

The PDF measures 17.9 cm × 11.6 cm and contains six line plots, a-f, in two rows and three columns.
The previous velocity maps have been removed. There are no subplot titles or embedded tables. Panel labels use bold 9 pt text without parentheses.
Other text uses 8 pt; compact annotations and the legend in the error-bound panel use 7.1 pt.

Colors and styling follow the supplied dense_sparse_consistency.pdf and its plotting notebook:
red #CE5044, blue #7AA7D3, green #86B6A4, and purple #A486BE, with corresponding pale marker fills.
A coordinated ochre #D6A55B extends the palette to the fifth grid.
Curves use thin solid lines, and marker edges match the line colors. All legends have white backgrounds and light-gray rectangular frames.
The panel-a legend uses the actual curves with markevery, without line-only proxy artists.

Run all cells in plot_mps_sixpanel.ipynb or execute the Python script with the same base name.
The reference styling requires NumPy, pandas, Matplotlib, and LaTeX/dvipng with sansmathfonts.
All inputs reside in data/; paths to the original project and reference figure are unnecessary.

Outputs include the PDF, 600 dpi PNG, figure_caption.tex, figure_include.tex,
two separate LaTeX tables, and captions.md. Caption references use panels a-f.
The two tables and all numerical plotting inputs match the previous version; this plotting revision did not regenerate fields or run MPS calculations.
Legacy slices and metadata remain in data/ for provenance and are not plotted in the new figure.

The circuit table includes the 30-qubit row, obtained by parameter compilation and transpilation without field or MPS simulation.
layout_validation.json records font sizes, legend correspondence, physical dimensions, and boundary checks.
validation.json records final-PDF and data-integrity checks.
