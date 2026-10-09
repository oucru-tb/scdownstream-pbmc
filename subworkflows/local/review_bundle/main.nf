include { ADATA_SETCELLTYPE  } from '../../../modules/local/adata/setcelltype'
include { ADATA_REVIEWBUNDLE } from '../../../modules/local/adata/reviewbundle'

workflow REVIEW_BUNDLE {
    take:
    ch_h5ad        // channel: [ meta, h5ad ]
    cell_type_col  //   value: string
    source_cols    //   value: string (space-separated, priority-ordered)
    group_cols     //   value: string (space-separated)

    main:
    ch_versions = channel.empty()

    ADATA_SETCELLTYPE (
        ch_h5ad,
        source_cols,
        cell_type_col,
    )
    ch_versions = ch_versions.mix(ADATA_SETCELLTYPE.out.versions)

    ADATA_REVIEWBUNDLE (
        ADATA_SETCELLTYPE.out.h5ad,
        cell_type_col,
        group_cols,
    )
    ch_versions = ch_versions.mix(ADATA_REVIEWBUNDLE.out.versions)

    emit:
    h5ad            = ADATA_SETCELLTYPE.out.h5ad  // channel: [ meta, h5ad ]
    review_h5ad     = ADATA_REVIEWBUNDLE.out.h5ad // channel: [ meta, h5ad ]
    annotation_csvs = ADATA_REVIEWBUNDLE.out.annotation_csvs
    plots           = ADATA_REVIEWBUNDLE.out.plots
    versions        = ch_versions
}
