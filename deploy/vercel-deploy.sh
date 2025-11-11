#!/bin/bash

# CriOS Discovery Engine - Vercel Deployment Script
# Quick deployment helper for Vercel

set -e

echo "🔷 CriOS Discovery Engine - Vercel Deployment Helper"
echo "=================================================="
echo ""

# Check if Vercel CLI is installed
if ! command -v vercel &> /dev/null; then
    echo "❌ Vercel CLI not found. Installing..."
    echo ""
    npm install -g vercel
    echo ""
fi

echo "✅ Vercel CLI found"
echo ""

# Check if logged in
if ! vercel whoami &> /dev/null; then
    echo "🔐 Please login to Vercel..."
    vercel login
fi

echo "✅ Logged in to Vercel"
echo ""

# Deployment choice
echo "Choose deployment option:"
echo ""
echo "  1) Hybrid (Vercel Frontend + Railway Backend) - RECOMMENDED"
echo "  2) Pure Vercel (Frontend + API Routes) - Limited chemistry features"
echo ""
read -p "Enter choice (1 or 2): " deploy_choice

if [ "$deploy_choice" = "1" ]; then
    echo ""
    echo "🏗️  Hybrid Deployment Selected"
    echo ""
    echo "This will deploy:"
    echo "  ✓ Frontend to Vercel (Next.js UI)"
    echo "  ✓ Backend to Railway (FastAPI + RDKit)"
    echo ""

    read -p "Have you already deployed the backend to Railway? (y/n): " backend_deployed

    if [ "$backend_deployed" != "y" ]; then
        echo ""
        echo "⚠️  Please deploy the backend first:"
        echo "   ./deploy/railway-deploy.sh"
        echo ""
        echo "Then run this script again."
        exit 1
    fi

    read -p "Enter your Railway backend URL (e.g., https://backend-production-xxxx.up.railway.app): " backend_url

    if [ -z "$backend_url" ]; then
        echo "❌ Backend URL is required for hybrid deployment"
        exit 1
    fi

    echo ""
    echo "📦 Deploying frontend to Vercel..."
    echo ""

    cd platform/frontend

    # Set environment variables
    vercel env add NEXT_PUBLIC_API_URL production <<< "$backend_url"
    vercel env add NEXT_PUBLIC_WS_URL production <<< "${backend_url/https/wss}"

    # Deploy
    vercel --prod

    echo ""
    echo "✅ Frontend deployed to Vercel!"
    echo ""
    echo "🔧 Don't forget to update your Railway backend CORS settings:"
    echo "   CORS_ORIGINS=https://your-project.vercel.app"
    echo ""

elif [ "$deploy_choice" = "2" ]; then
    echo ""
    echo "🏗️  Pure Vercel Deployment Selected"
    echo ""
    echo "⚠️  Note: This deployment has limited chemistry features"
    echo "   RDKit is too large for Vercel serverless functions"
    echo ""

    read -p "Continue with pure Vercel deployment? (y/n): " continue_pure

    if [ "$continue_pure" != "y" ]; then
        echo "Deployment cancelled."
        exit 0
    fi

    echo ""
    echo "📋 You'll need to set up:"
    echo "  1. Vercel Postgres or external database (Neon, Supabase)"
    echo "  2. Upstash Redis for caching"
    echo "  3. Anthropic API key"
    echo ""

    read -p "Have you set up external services? (y/n): " services_ready

    if [ "$services_ready" != "y" ]; then
        echo ""
        echo "⚠️  Please set up services first:"
        echo "   - Vercel Postgres: https://vercel.com/docs/storage/vercel-postgres"
        echo "   - Upstash Redis: https://console.upstash.com/"
        echo ""
        echo "Then run this script again."
        exit 1
    fi

    cd platform/frontend

    # Link or create project
    vercel

    echo ""
    echo "📝 Now add environment variables in Vercel Dashboard:"
    echo "   https://vercel.com/dashboard → Settings → Environment Variables"
    echo ""
    echo "Required variables:"
    echo "  - NEXT_PUBLIC_API_URL"
    echo "  - REDIS_URL"
    echo "  - POSTGRES_URL"
    echo "  - ANTHROPIC_API_KEY"
    echo ""

    read -p "Press Enter after adding environment variables..."

    # Deploy to production
    vercel --prod

    echo ""
    echo "✅ Deployed to Vercel!"
    echo ""

else
    echo "❌ Invalid choice"
    exit 1
fi

echo ""
echo "🎉 Deployment complete!"
echo ""
echo "📊 View deployment:"
echo "   vercel open"
echo ""
echo "📝 View logs:"
echo "   vercel logs --follow"
echo ""
echo "📚 For more details, see VERCEL_DEPLOYMENT.md"
echo ""
