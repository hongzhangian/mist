{
  'variables': {
    'mist_product_name': 'Mist',
    'mist_version': '<!(python build/version.py -f VERSION -t "@MAJOR@.@MINOR@.@BUILD@.@PATCH@")',
  },
  'targets': [
    {
      'target_name': 'mist_runtime',
      'type': 'shared_library',
      'defines': ['MIST_VERSION="<(mist_version)"'],
      'dependencies': [
      ],
      'include_dirs': [
        '.',
      ],
      'sources': [
      ],
    },
  ], # targets
}
