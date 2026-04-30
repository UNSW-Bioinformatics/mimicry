from shiny import App, render, ui, reactive
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import shinywidgets as sw

# Import data from shared.py
from shared import df, pairs
from covfiles import covlist, pcov, world, regions, countries, getPoptype, matchData, totalCov, covstring

# Import data for differential abundance
from differential_abundance_data import inhouse_hcc_data as hcc
from differential_abundance_data import bmi_dataframe as bmi

# Define the color palette
color_palette = ['#f2c1a1', '#8cc8e3', "#075e7a", "#262a54", '#d6a9e1', '#9b7bc6', '#e5a3c9', '#6e4b75', '#a5e4c9', '#3e5940']

# Hard-coded data for some charts
STARTING_SAMPLES = {
    'Bladder': 677, 'Brain': 243, 'Breast': 1769, 'Cervix': 425, 'Colorectal': 698,
    'Head & Neck': 780, 'Kidney': 852, 'Liver': 537, 'Lung': 1715, 'Ovary': 621,
    'Pancreas': 238, 'Prostate': 874, 'Skin': 753, 'Stomach': 589, 'Thyroid': 598, 'Uterus': 613
}

DATABASE_COUNTS = {'TSNAdb': 11982}

ALL_CANCERS = [
    'Bladder', 'Brain', 'Breast', 'Cervix', 'Colorectal', 'Head & Neck',
    'Kidney', 'Liver', 'Lung', 'Ovary', 'Pancreas', 'Prostate',
    'Skin', 'Stomach', 'Thyroid', 'Uterus'
]

# Function to get unique values for Cancer
def get_unique_values(column_name):
    return df[column_name].unique().tolist()

# Define fixed colors for each cancer type
def get_color_mapping(cancer_types):
    color_mapping = {}
    for i, cancer in enumerate(cancer_types):
        color_mapping[cancer] = color_palette[i % len(color_palette)]
    return color_mapping

# UI definition
app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.panel_conditional(
            "input.tabs == 'database_overview', 'mutation_analysis', 'body_location', 'microbiome'",
            ui.input_checkbox_group(
                "cancerType",
                "Select Cancer Types",
                choices=ALL_CANCERS,
                selected=ALL_CANCERS
            ),
            ui.div(
                ui.input_action_button("select_all", "Select all"),
                ui.input_action_button("select_none", "Deselect all"),
                class_="d-flex gap-2 flex-wrap"
            ),
            ui.input_checkbox_group(
                "mutationType",
                "Select Mutation Types",
                choices=["SNV", "Fusion", "INDEL"],
                selected=["SNV", "Fusion", "INDEL"]
            ),
            ui.input_switch("pidentFilter", "Show 100.00", value=False),
            ui.input_action_button("reset", "Reset filter")
        ),
        ui.panel_conditional(
            "input.tabs == 'diff_abundance'",
            ui.input_select("comparison_type", "Select a comparison", choices=["Liver disease vs Healthy", "BMI comparison"]),
            ui.panel_conditional(
                "input.comparison_type == 'Liver disease vs Healthy'",
                ui.input_select("sample", "Select sample", choices=["Oral", "Stool"]),
                ui.input_checkbox_group(
                    "disease_group",
                    "Disease type",
                    {
                        "NAFLD-cirrhosis": "NAFLD-cirrhosis",
                        "NAFLD-HCC": "NAFLD-HCC",
                        "non-NAFLD-HCC": "non-NAFLD-HCC",
                    },
                    selected="NAFLD-cirrhosis"
                ),
                ui.input_slider(id="inhouse_pvalue_slider", label="P-value threshold", min=0.0001, max=1, value=0.0500, step=0.0001),
                ui.input_slider(id="inhouse_fc_slider", label="Fold change threshold", min=1, max=50, value=1.5, step=0.1)
            ),
            ui.panel_conditional(
                "input.comparison_type == 'BMI comparison'",
                ui.input_select("cancer_group", "Select location", choices=["Liver", "Breast", "Kidney"]),
                ui.input_checkbox_group(
                    "bmi_group",
                    "BMI group",
                    {
                        "Overweight": "Overweight vs healthy weight",
                        "Obesity": "Obesity vs healthy weight",
                    },
                    selected="Overweight"
                ),
                ui.input_slider(id="bmi_pvalue_slider", label="P-value threshold", min=0.0001, max=1, value=0.0500, step=0.0001),
                ui.input_slider(id="bmi_fc_slider", label="Fold change threshold", min=1, max=50, value=1.5, step=0.1)
            ),
        ),
        ui.panel_conditional(
            "input.tabs == 'pop_coverage'",
            ui.input_select('cancer', 'Cancer', choices=['Select a body site'] + sorted(list(set(covlist['Cancer'])))),
            ui.input_select('ncombs', 'Number of Epitopes', choices=list(range(0,31))),
            ui.navset_card_underline(
                ui.nav_panel("Worldwide"),
                ui.nav_panel("by Region", ui.input_select('regi', '', choices=['Select a region'] + list(regions['Region']))),
                ui.nav_panel("by Country", ui.input_select('country', '', choices=['Select a country'] + list(countries['Country']))),
                title="Population", id="population",
            ),
        ),
        id='mainsidebar',
    ),
    ui.card(
        ui.navset_tab(
            ui.nav_panel("Home",                
                ui.output_image("figure", height="100px"),
                ui.card(                    
                    ui.markdown("""
                    ## MIMICRY: A database of Microbial Immunogenic Mimicry in Cancer
                    ### Overview
                    * Recent studies suggest a connection between the human microbiome and cancer, particularly through immune modulation.
                    * Here, we present a database of eptitopes that are shared between common cancers and the human microbiome from HGMD.
                    * Tumor antigens have been compared to all peptides in all human microbiota generating a list of shared epitopes.  
                    * Our resource aims to facilitate the exploration of shared epitopes in cancer research, potentially leading to breakthroughs in understanding tumor immunity and the development of novel immunotherapeutic strategies.
                    """)
                ),  
                value = "home"
            ),
            ui.nav_panel("Database Overview",
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Distribution of Variants",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Shows the proportion of different variant types (SNV, Fusion, INDEL) in result dataset.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("variant_pie_chart")
                ),
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Starting Samples vs Full Matches",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Compares the initial number of samples per cancer type to the results returned per cancer type.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("samples_vs_matches_bar")
                ),
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Data per Database",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Shows the total number of entries in each database source. Currently displays data from TSNAdb.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("database_bar")
                ),
                value = "database_overview"
            ),
            ui.nav_panel("Mutation Analysis",
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Number of Mutations per Cancer",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Total count of mutations detected in each cancer type. Higher counts indicate more mutation events.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("mutation_types_per_cancer_bar")
                ),
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Samples per Cancer / Variant",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Stacked percentage bar chart showing the proportion of each variant type within each cancer. Bars are normalised to 100%, with hover showing actual counts.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("samples_per_cancer_bar")
                ),
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Unique Queries per Cancer / Variant",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Shows the number of unique queries with each variant type per cancer.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("individuals_per_cancer_bar")
                ),
                value = "mutation_analysis"
            ),
            ui.nav_panel("Body Location Analysis",
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Counts for Each Body Location",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Displays the number of microbial mutations found at different body locations (oral, gut, skin, etc.).",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("body_bar_plot")
                ),
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Mutations per Body Location",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Bar chart summarizing the total number of mutations detected at each body location (oral, gut, skin, etc.) across all selected filters.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("mutation_body_location_bar")
                ),
                value = "body_location"
            ),
            ui.nav_panel("Microbiome Characterization",
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Counts for Each Class",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Shows the distribution of mutation classes across selected cancer types. Each bar represents a different cancer type, grouped by their classification.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("class_bar_plot")
                ),
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Epitope Binding Affinity Heatmap",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Heatmap of the top 20 epitopes (by max Deep bind score) across different microbe locations. Lighter colors indicate stronger predicted binding affinity (0-1 scale).",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("epitope_heatmap")
                ),
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Mutation Load per Tissue Type",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Scatter plot showing the number of neoantigens per cancer type, with color indicating average Deep immunogenicity score.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("mutation_load_scatter")
                ),
                ui.card(
                    ui.card_header(
                        ui.div(
                            "Binding Affinity by Mutation Type",
                            ui.popover(
                                ui.span("ℹ️", style="cursor: pointer; margin-left: 10px;"),
                                "Violin plots showing the distribution of Net4 binding affinity (nM) for different mutation types. Lower values indicate stronger binding.",
                                placement="right"
                            ),
                            style="display: flex; align-items: center;"
                        )
                    ),
                    sw.output_widget("binding_affinity_density")
                ),
                value = "microbiome"
            ),
            ui.nav_panel("Differential Abundance",
                ui.panel_conditional(
                    "input.comparison_type == 'Liver disease vs Healthy'",
                    ui.output_text('show_hcc_description'),
                    sw.output_widget("inhouse_liver_volcano_plot"),
                    ui.output_data_frame("show_df")
                ),
                ui.panel_conditional(
                    "input.comparison_type == 'BMI comparison'",
                    ui.output_text('show_bmi_description'),
                    sw.output_widget("bmi_volcano_plot"),
                    ui.output_data_frame("show_df_bmi"),
                ),
                value = "diff_abundance"
            ),
            ui.nav_panel("Genome Locations", ui.markdown("This is the third page.")),
            ui.nav_panel("Epitope Coverage",
                ui.card(
                    ui.output_text('cov_summary'),
                    ui.h5(ui.output_text('coverage')),
                    ui.output_data_frame('values'),
                ),
                ui.card(ui.output_plot('cov_bar_plot', height='500px')),
                value = "pop_coverage"
            ),
            id='tabs',
        ),
    )
)

def server(input, output, session):

    # --- Helper function to normalize SNP and SNV entries
    def normalize_variant(variant_series):
        return variant_series.replace('SNP', 'SNV')

    # --- Reactive function to get filtered dataframe
    @reactive.Calc
    def get_filtered_df():
        selected_cancers = input.cancerType()
        selected_mutations = input.mutationType()
        show_100_only = input.pidentFilter()
        
        filtered = df.copy()
        filtered['Variant'] = normalize_variant(filtered['Variant'])
        
        # Filter by cancer type
        if len(selected_cancers) > 0:
            filtered = filtered[filtered['Cancer'].isin(selected_cancers)]
        
        # Filter by mutation type
        if len(selected_mutations) > 0:
            filtered = filtered[filtered['Variant'].isin(selected_mutations)]
        
        # Filter by PIdentity and Coverage if toggle is on (> 99.9 for both)
        if show_100_only:
            filtered['PIdentity_num'] = pd.to_numeric(filtered['PIdentity'], errors='coerce')
            filtered['Coverage_num'] = pd.to_numeric(filtered['Coverage'], errors='coerce')
            
            # Filter for values > 99.9
            filtered = filtered[
                (filtered['PIdentity_num'] > 99.9) & 
                (filtered['Coverage_num'] > 99.9)
            ]
        
        return filtered

    # --- Buttons to (un)select cancers in the sidebar
    @reactive.Effect
    @reactive.event(input.select_none)
    def _untick_all():
        ui.update_checkbox_group("cancerType", selected=[])

    @reactive.Effect
    @reactive.event(input.select_all)
    def _select_all():
        ui.update_checkbox_group("cancerType", selected=ALL_CANCERS)

    # Reset button
    @reactive.Effect
    @reactive.event(input.reset)
    def _reset():
        ui.update_checkbox_group("cancerType", selected=ALL_CANCERS)
        ui.update_checkbox_group("mutationType", selected=["SNV", "Fusion", "INDEL"])
        ui.update_switch("pidentFilter", value=False)

    # 1. Class bar plot
    @output
    @sw.render_widget
    def class_bar_plot():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Counts for Each Class (no data matching filters)")
        class_counts = filtered_df.pivot_table(index='Class', columns='Cancer', aggfunc='size', fill_value=0)
        color_mapping = get_color_mapping(get_unique_values('Cancer'))
        fig = go.Figure()
        for cancer in class_counts.columns:
            fig.add_trace(go.Bar(
                name=cancer, x=class_counts.index, y=class_counts[cancer],
                marker_color=color_mapping.get(cancer, '#333333'),
                hovertemplate='<b>%{x}</b><br>Count: %{y}<extra></extra>'
            ))
        fig.update_layout(title="Counts for Each Class", xaxis_title="Class", yaxis_title="Count", barmode='group', hovermode='closest')
        return fig

    # 2. Body location bar plot
    @output
    @sw.render_widget
    def body_bar_plot():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Counts for Each Body Location (no data matching filters)")
        body_counts = filtered_df['Microbe location'].value_counts()
        fig = go.Figure(data=[go.Bar(x=body_counts.index, y=body_counts.values, marker_color=color_palette[3],
                                     hovertemplate='<b>%{x}</b><br>Count: %{y}<extra></extra>')])
        fig.update_layout(title="Counts for Each Body Location", xaxis_title="Body Location", yaxis_title="Count", hovermode='closest')
        return fig

    # 4. Distribution of Variants pie chart
    @output
    @sw.render_widget
    def variant_pie_chart():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Distribution of Variants (no data matching filters)")
        variant_counts = filtered_df['Variant'].value_counts()
        fig = go.Figure(data=[go.Pie(labels=variant_counts.index, values=variant_counts.values,
                                     marker=dict(colors=color_palette),
                                     hovertemplate='<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>')])
        fig.update_layout(title="Distribution of Variants (SNV/Indel/Fusion)")
        return fig

    # 5. Starting samples vs full matches
    @output
    @sw.render_widget
    def samples_vs_matches_bar():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Starting Samples vs Full Matches by Cancer Type (no data matching filters)")
        filtered_df_copy = filtered_df.copy()
        filtered_df_copy['PIdentity_num'] = pd.to_numeric(filtered_df_copy['PIdentity'], errors='coerce')
        full_matches = filtered_df_copy[filtered_df_copy['PIdentity_num'] >= 99.999].groupby('Cancer').size()
        cancers = sorted(list(set(filtered_df['Cancer'])))
        starting = [STARTING_SAMPLES.get(c, 0) for c in cancers]
        matches = [full_matches.get(c, 0) for c in cancers]
        fig = go.Figure(data=[
            go.Bar(name='Starting Samples', x=cancers, y=starting, marker_color=color_palette[0],
                   hovertemplate='<b>%{x}</b><br>Starting: %{y}<extra></extra>'),
            go.Bar(name='Full Matches', x=cancers, y=matches, marker_color=color_palette[1],
                   hovertemplate='<b>%{x}</b><br>Matches: %{y}<extra></extra>')
        ])
        fig.update_layout(title="Starting Samples vs Full Matches by Cancer Type", xaxis_title="Cancer Type", yaxis_title="Count", barmode='group')
        return fig

    # 6. Number of mutations per cancer
    @output
    @sw.render_widget
    def mutation_types_per_cancer_bar():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Number of Mutations per Cancer (no data matching filters)")
        mutation_counts = filtered_df.groupby('Cancer').size().reset_index(name='count')
        fig = go.Figure(data=[go.Bar(x=mutation_counts['Cancer'], y=mutation_counts['count'],
                                     marker_color=color_palette[0],
                                     hovertemplate='<b>%{x}</b><br>Count: %{y}<extra></extra>')])
        fig.update_layout(title="Number of Mutations per Cancer", xaxis_title="Cancer", yaxis_title="Count", showlegend=False)
        return fig

    # 7. Samples per cancer/variant
    @output
    @sw.render_widget
    def samples_per_cancer_bar():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Number of Samples per Cancer / Variant (no data matching filters)")
        
        # Group by Cancer and Variant to get counts
        samples_counts = filtered_df.groupby(['Cancer', 'Variant']).size().reset_index(name='count')
        
        # Calculate total per cancer for percentage calculation
        cancer_totals = samples_counts.groupby('Cancer')['count'].sum().reset_index(name='total')
        samples_counts = samples_counts.merge(cancer_totals, on='Cancer')
        samples_counts['percentage'] = (samples_counts['count'] / samples_counts['total']) * 100
        
        fig = go.Figure()
        
        # Get unique variants for consistent coloring
        variants = samples_counts['Variant'].unique()
        variant_colors = {variant: color_palette[i % len(color_palette)] for i, variant in enumerate(variants)}
        
        # Create stacked bar chart with percentages
        for variant in variants:
            variant_data = samples_counts[samples_counts['Variant'] == variant]
            fig.add_trace(go.Bar(
                name=variant,
                x=variant_data['Cancer'],
                y=variant_data['percentage'],
                marker_color=variant_colors[variant],
                customdata=variant_data['count'],
                hovertemplate='<b>%{x}</b><br>%{fullData.name}<br>Count: %{customdata}<br>Percentage: %{y:.1f}%<extra></extra>'
            ))
        
        fig.update_layout(
            title="Number of Samples per Cancer / Variant",
            xaxis_title="Cancer",
            yaxis_title="Percentage (%)",
            barmode='stack',
            yaxis=dict(range=[0, 100])
        )
        return fig

    # 8. Individuals per cancer/variant
    @output
    @sw.render_widget
    def individuals_per_cancer_bar():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Number of Queries per Cancer / Variant (no data matching filters)")
        individuals_counts = filtered_df.groupby(['Cancer', 'Variant'])['Query'].nunique().reset_index(name='count')
        fig = px.bar(individuals_counts, x='Cancer', y='count', color='Variant',
                     title="Number of Queries per Cancer / Variant",
                     color_discrete_sequence=color_palette, labels={'count': 'Unique Individuals'})
        fig.update_traces(hovertemplate='<b>%{x}</b><br>%{fullData.name}: %{y}<extra></extra>')
        return fig

    # 9. Heatmap of epitopes with binding affinity (TOP 50, NaN-safe everywhere)
    @output
    @sw.render_widget
    def epitope_heatmap():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Heatmap (no data matching filters)")

        col = 'Deep bind'
        filtered_df[col] = (
            filtered_df[col].astype(str).str.strip().str.replace(',', '.', regex=False)
        )
        filtered_df[col] = pd.to_numeric(filtered_df[col], errors='coerce')
        scored = filtered_df.copy()
        scored = scored[~scored[col].isin([np.inf, -np.inf])]
        valid_scored = scored[(scored[col].notna()) & (scored[col] >= 0) & (scored[col] <= 1)]
        
        if valid_scored.empty:
            return go.Figure(layout_title_text="Heatmap (no valid Deep bind values)")

        # Top 20 queries by max Deep bind
        top_queries = (
            valid_scored.groupby('Query', as_index=False)[col]
                .max()
                .sort_values(col, ascending=False)
                .head(20)['Query']
        )
        top_df = scored[scored['Query'].isin(top_queries)]

        heatmap_df = top_df.pivot_table(
            index='Query',
            columns='Microbe location',
            values=col,
            aggfunc='max'
        )
        heatmap_df = heatmap_df.replace([np.inf, -np.inf], np.nan)
        
        # Convert to list of lists, replacing NaN and values outside range with 0 (purple)
        z_matrix = []
        hover_text = []
        for i, row in enumerate(heatmap_df.values):
            clean_row = []
            hover_row = []
            for j, val in enumerate(row):
                if pd.isna(val):
                    clean_row.append(0.0)
                    hover_row.append(f"Location: {heatmap_df.columns[j]}<br>Epitope: {heatmap_df.index[i]}<br>Deep bind: N/A")
                elif val < 0 or val > 1:
                    # Out of range values also appear as 0/purple
                    clean_row.append(0.0)
                    hover_row.append(f"Location: {heatmap_df.columns[j]}<br>Epitope: {heatmap_df.index[i]}<br>Deep bind: N/A (out of range)")
                else:
                    clean_row.append(float(val))
                    hover_row.append(f"Location: {heatmap_df.columns[j]}<br>Epitope: {heatmap_df.index[i]}<br>Deep bind: {val:.4f}")
            z_matrix.append(clean_row)
            hover_text.append(hover_row)

        def _labelize(vals):
            out = []
            for v in vals:
                if v is None or (isinstance(v, float) and pd.isna(v)):
                    out.append("Unknown")
                else:
                    out.append(str(v))
            return out

        x_vals = _labelize(heatmap_df.columns.tolist())
        y_vals = _labelize(heatmap_df.index.tolist())

        fig = go.Figure(data=go.Heatmap(
            z=z_matrix,
            x=x_vals,
            y=y_vals,
            colorscale='Viridis',
            zmin=0.0,
            zmax=1.0,
            hovertext=hover_text,
            hoverinfo='text',
            colorbar=dict(tickformat=".4f")
        ))
        fig.update_layout(
            title="Heatmap of Top 20 Epitopes (by max Deep bind)",
            xaxis_title="Microbe Location",
            yaxis_title="Epitope (Query)",
            height=700
        )
        return fig

    # 10. Mutation load per tissue type
    @output
    @sw.render_widget
    def mutation_load_scatter():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Mutation Load per Tissue Type (no data matching filters)")
        mutation_load = filtered_df.groupby('Cancer').agg({'Query': 'count', 'Deep imm': 'mean'}).reset_index()
        mutation_load.columns = ['Cancer', 'Neoantigens', 'Avg_Deep_imm']
        fig = go.Figure(data=go.Scatter(
            x=mutation_load['Cancer'], y=mutation_load['Neoantigens'], mode='markers',
            marker=dict(size=15, color=mutation_load['Avg_Deep_imm'], colorscale='RdYlGn', showscale=True,
                        colorbar=dict(title="Avg Deep Imm")),
            text=mutation_load['Cancer'],
            hovertemplate='<b>%{text}</b><br>Neoantigens: %{y}<br>Avg Deep Imm: %{marker.color:.3f}<extra></extra>'
        ))
        fig.update_layout(title="Mutation Load per Tissue Type", xaxis_title="Tissue Type", yaxis_title="Number of Neoantigens")
        return fig

    # 11. Database bar
    @output
    @sw.render_widget
    def database_bar():
        databases = list(DATABASE_COUNTS.keys())
        counts = list(DATABASE_COUNTS.values())
        fig = go.Figure(data=[go.Bar(x=databases, y=counts, marker_color=color_palette[2],
                                     hovertemplate='<b>%{x}</b><br>Count: %{y}<extra></extra>')])
        fig.update_layout(title="Data per Database", xaxis_title="Database", yaxis_title="Count")
        return fig

    # 12. Binding affinity density
    @output
    @sw.render_widget
    def binding_affinity_density():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Mutation Type vs Binding Affinity (Net4_aff) (no data matching filters)")
        plot_df = filtered_df[filtered_df['Net4 aff (nM)'].notna()]
        fig = go.Figure()
        for i, mutation in enumerate(plot_df['Mutation'].unique()):
            mutation_data = plot_df[plot_df['Mutation'] == mutation]['Net4 aff (nM)']
            fig.add_trace(go.Violin(
                y=mutation_data, name=mutation, box_visible=True, meanline_visible=True,
                marker_color=color_palette[i % len(color_palette)],
                hovertemplate='%{fullData.name}<br>Value: %{y:.2f}<extra></extra>'
            ))
        fig.update_layout(title="Mutation Type vs Binding Affinity (Net4_aff)",
                          yaxis_title="Net4 aff (nM)", xaxis_title="Mutation Type", showlegend=False)
        return fig

    # 13. Mutations per body location
    @output
    @sw.render_widget
    def mutation_body_location_bar():
        filtered_df = get_filtered_df()
        if filtered_df.empty:
            return go.Figure(layout_title_text="Mutations per Body Location (no data matching filters)")
        mutation_location = filtered_df.groupby('Microbe location').size().reset_index(name='count')
        fig = go.Figure(data=[go.Bar(x=mutation_location['Microbe location'], y=mutation_location['count'],
                                     marker_color=color_palette[2],
                                     hovertemplate='<b>%{x}</b><br>Count: %{y}<extra></extra>')])
        fig.update_layout(title="Mutations per Body Location", xaxis_title="Body Location", yaxis_title="Count", showlegend=False)
        return fig

    @render.text
    def show_hcc_description():
        return "This page shows differential abundance of microbial in healthy people vs patients with " + ", ".join(input.disease_group()) + " in either oral or stool. The positive fold changes indicate more abundance in patients compare to healthy people. The data frame under the plot shows the calculated data that is used to show the abundance."

    @output
    @sw.render_widget
    def inhouse_liver_volcano_plot():
        df_sample = hcc.df_hcc_oral if input.sample() == "Oral" else hcc.df_hcc_stool
        fc_threshold = np.log2(input.inhouse_fc_slider())
        pvalue_threshold = np.log10(input.inhouse_pvalue_slider()) * (-1)
        fig = go.Figure()
        disease_type = input.disease_group()
        if "NAFLD-cirrhosis" in disease_type:
            fig.add_trace(go.Scatter(x=df_sample['log2FC_OCI'], y=df_sample['negative_log_pval_oci'],
                                     mode='markers', name='NAFLD-cirrhosis', hovertext=list(df_sample.index)))
        if "NAFLD-HCC" in disease_type:
            fig.add_trace(go.Scatter(x=df_sample['log2FC_OLN'], y=df_sample['negative_log_pval_oln'],
                                     mode='markers', name='NAFLD-HCC', hovertext=list(df_sample.index)))
        if "non-NAFLD-HCC" in disease_type:
            fig.add_trace(go.Scatter(x=df_sample['log2FC_OLX'], y=df_sample['negative_log_pval_olx'],
                                     mode='markers', name='Non-NAFLD-HCC', hovertext=list(df_sample.index)))
        fig.add_hline(y=pvalue_threshold, line_width=1, line_dash="dash", line_color="green")
        fig.add_vline(x=fc_threshold, line_width=1, line_dash="dash", line_color="blue")
        fig.add_vline(x=(fc_threshold * (-1)), line_width=1, line_dash="dash", line_color="blue")
        fig.update_layout(xaxis_title="log2FC", yaxis_title="-log10(pvalue)")
        return fig

    @output
    @render.data_frame
    def show_df():
        if input.sample() == "Oral":
            df_sample = hcc.df_hcc_oral
        elif input.sample() == "Stool":
            df_sample = hcc.df_hcc_stool
        selected_columns = df_sample[['organism', 'mean OCO']]
        disease_type = input.disease_group()
        if "NAFLD-cirrhosis" in disease_type:
            selected_columns = pd.concat([selected_columns, df_sample[['mean OCI', 'pvalue OCI to control', 'ratio OCI to control', 'log2FC_OCI', 'negative_log_pval_oci']]], axis=1)
        if "NAFLD-HCC" in disease_type:
            selected_columns = pd.concat([selected_columns, df_sample[['mean OLX', 'pvalue OLX to control', 'ratio OLX to control', 'log2FC_OLX', 'negative_log_pval_olx']]], axis=1)
        if "non-NAFLD-HCC" in disease_type:
            selected_columns = pd.concat([selected_columns, df_sample[['mean OLN', 'pvalue OLN to control', 'ratio OLN to control', 'log2FC_OLN', 'negative_log_pval_oln']]], axis=1)
        return render.DataGrid(selected_columns)

    @render.text
    def show_bmi_description():
        return "This page shows differential abundance of microbial in healthy weight group vs " + ", ".join(input.bmi_group()) + "group. The organ location selection is from SNV data in TSNAdb. The positive fold changes indicate more abundance in treatment group compare to healthy weight group. The data frame under the plot shows the calculated data that is used to show the abundance."

    @output
    @sw.render_widget
    def bmi_volcano_plot():
        if input.cancer_group() == "Liver":
            df_sample = bmi.df_bmi_liver
        elif input.cancer_group() == "Breast":
            df_sample = bmi.df_bmi_breast
        elif input.cancer_group() == "Kidney":
            df_sample = bmi.df_bmi_kidney
        fc_threshold = np.log2(input.bmi_fc_slider())
        pvalue_threshold = np.log10(input.bmi_pvalue_slider()) * (-1)
        fig = go.Figure()
        bmi_type = input.bmi_group()
        if "Overweight" in bmi_type:
            fig.add_trace(go.Scatter(x=df_sample['log2FC_Overweight'], y=df_sample['negative_log_pval_overweight'],
                                     mode='markers', name='Overweight', hovertext=list(df_sample.index)))
        if "Obesity" in bmi_type:
            fig.add_trace(go.Scatter(x=df_sample['log2FC_Obesity'], y=df_sample['negative_log_pval_obesity'],
                                     mode='markers', name='Obesity', hovertext=list(df_sample.index)))
        fig.add_hline(y=pvalue_threshold, line_width=1, line_dash="dash", line_color="green")
        fig.add_vline(x=fc_threshold, line_width=1, line_dash="dash", line_color="blue")
        fig.add_vline(x=(fc_threshold * (-1)), line_width=1, line_dash="dash", line_color="blue")
        fig.update_layout(xaxis_title="log2FC", yaxis_title="-log10(pvalue)")
        return fig

    @output
    @render.data_frame
    def show_df_bmi():
        if input.cancer_group() == "Liver":
            df_sample = bmi.df_bmi_liver
        elif input.cancer_group() == "Breast":
            df_sample = bmi.df_bmi_breast
        elif input.cancer_group() == "Kidney":
            df_sample = bmi.df_bmi_kidney
        selected_columns = df_sample[['organism', 'mean healthy weight']]
        bmi_type = input.bmi_group()
        if "Overweight" in bmi_type:
            selected_columns = pd.concat([selected_columns, df_sample[['mean overweight', 'pvalue Overweight to healthy weight', 'ratio of Overweight to healthy weight', 'log2FC_Overweight', 'negative_log_pval_overweight']]], axis=1)
        if "Obesity" in bmi_type:
            selected_columns = pd.concat([selected_columns, df_sample[['mean obesity', 'pvalue Obesity to healthy weight', 'ratio of Obesity to healthy weight', 'log2FC_Obesity', 'negative_log_pval_obesity']]], axis=1)
        return render.DataGrid(selected_columns)

    @output
    @render.plot
    def cov_bar_plot():
        area, param = getPoptype(input.population(), input.regi(), input.country(), regions, countries)
        coverage = area['Total Coverage %'].tolist()
        cancer = area['Cancer'].tolist()
        if param != 'World':
            df_plot = pd.DataFrame({param : area[param], 'Coverage' : coverage, 'Cancer' : cancer})
            rotate = 90 if param == 'Country' else 70
            df_plot.set_index([param, 'Cancer'], inplace=True)
            ax = df_plot.unstack().plot.bar(y='Coverage', ylabel='Predicted Total Coverage (%)', rot=rotate, color=color_palette)
            ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        else:
            df_plot = pd.DataFrame({'Cancer' : cancer, 'Coverage' : coverage})
            df_plot.plot.bar(x='Cancer', ylabel='Predicted Total Coverage (%)', rot=70, color=color_palette[2], legend=None)

    @output
    @render.data_frame
    def values():
        data = []
        pop = {'Worldwide' : 'World', 'by Region' : str(input.regi()), 'by Country' : str(input.country())}
        area = pop[input.population()]
        for i in range(0, len(covlist)):
            if covlist['Area'][i] == area and covlist['Cancer'][i] == str(input.cancer()):
                return matchData(data, i, int(input.ncombs()))

    @render.text
    def coverage():
        region, country, num, cancer, population = input.regi(), input.country(), str(input.ncombs()), input.cancer(), input.population()
        pop = {'Worldwide' : 'World', 'by Region' : region, 'by Country' : country}
        area = pop[population]
        retstr = 'Worldwide' if area == 'World' else f'in {area}'
        if not area.startswith('Select a') and cancer != 'Select a body site':
            ctype = pairs[cancer] if cancer in pairs.keys() else cancer
            if num == '0':
                return
            for i in range(0, len(pcov)):
                if pcov['Cancer'][i] == cancer and pcov['Area'][i] == area:
                    return totalCov(pcov, cancer, area, ctype, retstr, num, i)

    @render.text
    def cov_summary():
        return covstring
    
    @render.image
    def figure():
        from pathlib import Path
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / "figure.png"), "width": "100px"}
        return img
    
    @render.image
    def workflow():
        from pathlib import Path
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / "workflow.png"), "width": "200px"}

        return img
    
app = App(app_ui, server)

if __name__ == "__main__":
    app.run()