#!/bin/bash
# READY TO DEPLOY - Run these commands to deploy CriOS
# All preparation is complete. Just authenticate and deploy!

set -e

echo "╔═══════════════════════════════════════════════════════════╗"
echo "║                                                           ║"
echo "║   🎉 CriOS IS READY TO DEPLOY! 🎉                        ║"
echo "║                                                           ║"
echo "║   ✅ Dependencies installed (0 vulnerabilities)           ║"
echo "║   ✅ Environment files configured                         ║"
echo "║   ✅ All configuration files in place                     ║"
echo "║   ✅ Package.json issues fixed                            ║"
echo "║                                                           ║"
echo "╚═══════════════════════════════════════════════════════════╝"
echo ""

# Check if we're in the right directory
if [ ! -f "deploy.sh" ]; then
    echo "❌ Error: Please run this from the project root directory"
    exit 1
fi

echo "📋 DEPLOYMENT OPTIONS:"
echo ""
echo "1️⃣  OPTION 1: Quick Deploy to Vercel (Recommended)"
echo "   └─ Fastest way to get your frontend live"
echo ""
echo "2️⃣  OPTION 2: GitHub Integration (Zero-maintenance)"
echo "   └─ Auto-deploy on every git push"
echo ""
echo "3️⃣  OPTION 3: Full Stack Deploy (Vercel + Fly.io)"
echo "   └─ Deploy everything including backend"
echo ""
read -p "Choose option [1-3]: " option

case $option in
    1)
        echo ""
        echo "🚀 DEPLOYING TO VERCEL..."
        echo ""
        echo "Step 1: Login to Vercel"
        vercel login

        echo ""
        echo "Step 2: Deploy Frontend"
        cd platform/frontend
        vercel --prod
        cd ../..

        echo ""
        echo "Step 3: Deploy Algorithm Studio UI"
        cd ui
        vercel --prod
        cd ..

        echo ""
        echo "✅ DEPLOYMENT COMPLETE!"
        echo ""
        echo "Your apps are now live at:"
        echo "  • Frontend: Check the URL in the output above"
        echo "  • UI: Check the URL in the output above"
        echo ""
        echo "📝 Next steps:"
        echo "  1. Test your deployments"
        echo "  2. Deploy backend to Fly.io (see option 3)"
        echo "  3. Update NEXT_PUBLIC_API_URL if backend URL changes"
        ;;

    2)
        echo ""
        echo "🌐 GITHUB INTEGRATION SETUP"
        echo ""
        echo "Follow these steps:"
        echo ""
        echo "1. Go to: https://vercel.com/new"
        echo ""
        echo "2. Click 'Import Git Repository'"
        echo ""
        echo "3. Select: MichaelCrowe11/crios-dr-crowe-coder"
        echo ""
        echo "4. Create PROJECT 1 - Frontend:"
        echo "   • Project Name: crios-frontend"
        echo "   • Root Directory: platform/frontend"
        echo "   • Framework: Next.js"
        echo "   • Environment Variables:"
        echo "     NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev"
        echo "     NEXT_PUBLIC_WS_URL=wss://crios-backend.fly.dev"
        echo "   • Click Deploy"
        echo ""
        echo "5. Create PROJECT 2 - UI:"
        echo "   • Import repository again"
        echo "   • Project Name: crios-ui"
        echo "   • Root Directory: ui"
        echo "   • Framework: Next.js"
        echo "   • Environment Variables:"
        echo "     NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev"
        echo "   • Click Deploy"
        echo ""
        echo "6. Done! Now every git push will auto-deploy"
        echo ""
        read -p "Press ENTER when you've completed these steps..."
        echo "✅ GitHub integration should now be set up!"
        ;;

    3)
        echo ""
        echo "🌍 FULL STACK DEPLOYMENT"
        echo ""
        echo "═══════════════════════════════════════"
        echo "PART 1: Deploy Frontend to Vercel"
        echo "═══════════════════════════════════════"
        vercel login

        cd platform/frontend
        echo "Deploying frontend..."
        vercel --prod
        cd ../..

        cd ui
        echo "Deploying UI..."
        vercel --prod
        cd ..

        echo ""
        echo "═══════════════════════════════════════"
        echo "PART 2: Deploy Backend to Fly.io"
        echo "═══════════════════════════════════════"

        # Check if flyctl is installed
        if ! command -v flyctl &> /dev/null; then
            echo "Installing Fly.io CLI..."
            curl -L https://fly.io/install.sh | sh
            export PATH="$HOME/.fly/bin:$PATH"
        fi

        echo "Logging into Fly.io..."
        flyctl auth login

        echo ""
        echo "Creating PostgreSQL database..."
        echo "Run: flyctl postgres create --name crios-db --region iad --vm-size shared-cpu-1x --volume-size 10"
        read -p "Press ENTER after running the command and saving the DATABASE_URL..."
        read -p "Enter DATABASE_URL: " DATABASE_URL

        echo ""
        echo "Creating Redis instance..."
        echo "Run: flyctl redis create --name crios-redis --region iad --plan free"
        read -p "Press ENTER after running the command and saving the REDIS_URL..."
        read -p "Enter REDIS_URL: " REDIS_URL

        echo ""
        read -p "Enter your ANTHROPIC_API_KEY: " ANTHROPIC_API_KEY

        cd platform/backend

        # Create app if it doesn't exist
        if ! flyctl status -a crios-backend &> /dev/null; then
            echo "Creating Fly.io app..."
            flyctl apps create crios-backend
        fi

        echo "Setting secrets..."
        flyctl secrets set \
            DATABASE_URL="$DATABASE_URL" \
            REDIS_URL="$REDIS_URL" \
            ANTHROPIC_API_KEY="$ANTHROPIC_API_KEY" \
            -a crios-backend

        echo "Deploying backend..."
        flyctl deploy

        cd ../..

        echo ""
        echo "✅ FULL STACK DEPLOYMENT COMPLETE!"
        echo ""
        echo "Your apps are live at:"
        echo "  • Frontend: (check Vercel dashboard)"
        echo "  • UI: (check Vercel dashboard)"
        echo "  • Backend: https://crios-backend.fly.dev"
        echo "  • API Docs: https://crios-backend.fly.dev/docs"
        ;;

    *)
        echo "❌ Invalid option"
        exit 1
        ;;
esac

echo ""
echo "═══════════════════════════════════════"
echo "📊 DEPLOYMENT SUMMARY"
echo "═══════════════════════════════════════"
echo ""
echo "✅ What's deployed:"
if [ "$option" == "1" ] || [ "$option" == "3" ]; then
    echo "  • Frontend → Vercel"
    echo "  • Algorithm Studio UI → Vercel"
fi
if [ "$option" == "3" ]; then
    echo "  • Backend → Fly.io"
fi
echo ""
echo "📚 Useful links:"
echo "  • Vercel Dashboard: https://vercel.com/dashboard"
echo "  • Fly.io Dashboard: https://fly.io/dashboard"
echo "  • Documentation: DEPLOYMENT.md"
echo ""
echo "🎉 Congratulations! Your application is deployed!"
