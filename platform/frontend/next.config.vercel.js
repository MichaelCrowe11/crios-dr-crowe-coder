/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,

  // Image configuration
  images: {
    domains: ['localhost', 'crios.ai', 'vercel.app'],
    unoptimized: process.env.NODE_ENV === 'development',
  },

  // Environment variables
  env: {
    NEXT_PUBLIC_API_URL: process.env.NEXT_PUBLIC_API_URL || 'http://localhost:3000',
    NEXT_PUBLIC_WS_URL: process.env.NEXT_PUBLIC_WS_URL || 'ws://localhost:3000',
  },

  // Webpack configuration
  webpack: (config, { isServer }) => {
    // Handle canvas for server-side rendering
    config.externals = [...(config.externals || []), { canvas: 'canvas' }];

    // Ignore .py files in webpack bundling
    config.module.rules.push({
      test: /\.py$/,
      loader: 'ignore-loader',
    });

    return config;
  },

  // Experimental features for better performance
  experimental: {
    // Enable server actions (for Next.js 13+)
    serverActions: true,
  },

  // Headers for API routes
  async headers() {
    return [
      {
        source: '/api/:path*',
        headers: [
          { key: 'Access-Control-Allow-Credentials', value: 'true' },
          { key: 'Access-Control-Allow-Origin', value: '*' },
          { key: 'Access-Control-Allow-Methods', value: 'GET,DELETE,PATCH,POST,PUT' },
          { key: 'Access-Control-Allow-Headers', value: 'X-CSRF-Token, X-Requested-With, Accept, Accept-Version, Content-Length, Content-MD5, Content-Type, Date, X-Api-Version, Authorization' },
        ],
      },
    ];
  },

  // Rewrites for API routes
  async rewrites() {
    return [
      {
        source: '/health',
        destination: '/api/health',
      },
      {
        source: '/docs',
        destination: '/api/docs',
      },
    ];
  },

  // Redirects
  async redirects() {
    return [
      {
        source: '/api/docs',
        destination: '/api-docs',
        permanent: false,
      },
    ];
  },
}

module.exports = nextConfig
