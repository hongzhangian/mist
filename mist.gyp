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
        'math/math.gypi:mist_math'
      ],
      'include_dirs': [
        '.',
      ],
      'sources': [
        'math/vector2.cc',
        'math/vector2.h',
        'math/vector3.cc',
        'math/vector3.h',
        'math/vector4.cc',
        'math/vector4.h',
      ],
    },
  ], # targets
}
